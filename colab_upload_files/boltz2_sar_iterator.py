#!/usr/bin/env python3
"""
Boltz-2 SAR Iterator Tool

An iterative tool that runs Boltz-2 cofolding simulations for compounds with SAR data,
correlating predicted affinities with experimental activity values.
"""

import os
import json
import yaml
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import logging
from datetime import datetime
from scipy import stats
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

# Import PDB to CIF conversion utilities
try:
    import gemmi
    GEMMI_AVAILABLE = True
except ImportError:
    GEMMI_AVAILABLE = False
    logging.warning("Gemmi not available. PDB to CIF conversion will not work.")


def convert_pdb_to_cif(pdb_path: str, output_dir: Optional[str] = None) -> str:
    """
    Convert PDB file to mmCIF format using Gemmi.

    Args:
        pdb_path: Path to input PDB file
        output_dir: Optional output directory for CIF file

    Returns:
        Path to generated CIF file

    Raises:
        ImportError: If gemmi is not available
        FileNotFoundError: If PDB file doesn't exist
    """
    if not GEMMI_AVAILABLE:
        raise ImportError(
            "Gemmi is required for PDB to CIF conversion. "
            "Install with: pip install gemmi"
        )

    import re

    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    # Determine output path
    if output_dir:
        out_dir = Path(output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        cif_path = out_dir / f"{pdb_path.stem}.cif"
    else:
        cif_path = pdb_path.with_suffix('.cif')

    # Read structure
    st = gemmi.read_structure(str(pdb_path))

    # Chimera-like setup (from pdb2cifgemmi.py)
    def _sanitize_chain_names(st: gemmi.Structure):
        for model in st:
            for ch in model:
                name = (ch.name or "").strip()
                if not name or name in {".", "?"}:
                    ch.name = f"X{model.index}-{ch.index}"
                else:
                    ch.name = re.sub(r"[^A-Za-z0-9_-]", "_", name)

    def _populate_entity_sequences_from_model(st: gemmi.Structure):
        if len(st) == 0:
            return
        model0 = st[0]

        subchain_by_id = {}
        for ch in model0:
            for sc in ch.subchains():
                if hasattr(sc, "get_polymer"):
                    poly = sc.get_polymer()
                    sid = sc.subchain_id() if hasattr(sc, "subchain_id") else None
                else:
                    poly = sc
                    sid = sc.subchain_id() if hasattr(sc, "subchain_id") else None

                if sid is None:
                    try:
                        chain_name = ch.name
                    except Exception:
                        chain_name = "?"
                    sid = f"{chain_name}:{poly.first().seqid}" if len(poly) else f"{chain_name}:empty"

                if len(poly) == 0:
                    continue
                subchain_by_id[str(sid)] = poly

        for ent in st.entities:
            if ent.entity_type != gemmi.EntityType.Polymer:
                continue

            candidate_seqs = []
            for sid in ent.subchains:
                poly = subchain_by_id.get(str(sid))
                if poly is None or len(poly) == 0:
                    continue
                seq = [res.name for res in poly]
                candidate_seqs.append(seq)

            if not candidate_seqs:
                continue

            best = max(candidate_seqs, key=len)
            ent.full_sequence = best

    # Apply Chimera-like transformations
    _sanitize_chain_names(st)
    st.assign_subchains(force=True, fail_if_unknown=False)
    st.setup_entities()
    st.ensure_entities()
    st.add_entity_types(overwrite=True)
    st.add_entity_ids(overwrite=True)
    st.deduplicate_entities()
    _populate_entity_sequences_from_model(st)
    st.assign_label_seq_id(force=True)

    # Write mmCIF
    doc = st.make_mmcif_document()
    cif_path.write_text(doc.as_string())

    logging.info(f"Converted PDB to CIF: {cif_path}")
    return str(cif_path)


@dataclass
class SARData:
    """Container for Structure-Activity Relationship data."""
    compound_id: str
    smiles: str
    activity: float
    predicted_affinity: Optional[float] = None  # Primary metric for backward compatibility
    predicted_metrics: Optional[Dict[str, float]] = None  # All available metrics


class Boltz2SARIterator:
    """
    Iterative tool for running Boltz-2 cofolding simulations with SAR data.

    This tool:
    1. Loads SAR data from CSV (SMILES and Activity columns)
    2. Runs Boltz-2 cofolding simulations for each compound
    3. Extracts predicted affinities from Boltz-2 outputs
    4. Calculates correlation (R²) between predicted and experimental values
    5. Iterates until R² threshold or max iterations is reached
    """

    def __init__(
        self,
        protein_sequence: str,
        csv_path: str,
        output_dir: str = "./boltz2_sar_output",
        target_r2: Optional[float] = None,
        target_spearman: Optional[float] = None,
        target_roc_auc: Optional[float] = None,
        max_iterations: int = 10,
        use_msa_server: bool = True,
        msa_path: Optional[str] = None,
        protein_chain_id: str = "A",
        ligand_chain_id: str = "L",
        pocket_residues: Optional[List[Tuple[str, int]]] = None,
        contact_residues: Optional[List[Tuple[Tuple[str, int], Tuple[str, int], float]]] = None,
        template_file: Optional[str] = None,
        template_force: bool = True,
        template_threshold: int = 2,
        roc_threshold: Optional[float] = None,
        log_level: str = "INFO"
    ):
        """
        Initialize the Boltz-2 SAR Iterator.

        Args:
            protein_sequence: Protein amino acid sequence
            csv_path: Path to CSV file with SMILES and Activity columns
            output_dir: Directory for outputs and intermediate files
            target_r2: Target R² threshold (converges if any metric's R² >= this). Default: 0.6
            target_spearman: Target Spearman R threshold (converges if any metric >= this). Default: 0.7
            target_roc_auc: Target ROC-AUC threshold (converges if any metric >= this). Default: 0.8
            max_iterations: Maximum number of iterations
            use_msa_server: Whether to use MSA server for alignments
            msa_path: Optional path to pre-computed MSA file (.a3m)
            protein_chain_id: Chain ID for protein (default: 'A')
            ligand_chain_id: Chain ID for ligand (default: 'L')
            pocket_residues: List of (chain, residue) tuples defining binding pocket
            contact_residues: List of ((chain1, res1), (chain2, res2), max_dist) for contacts
            template_file: Path to template CIF file
            template_force: Force template usage (default: True)
            template_threshold: Template threshold (default: 2)
            roc_threshold: Optional threshold for ROC-AUC calculation (e.g., 10000 for IC50 < 10μM = binder)
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        self.protein_sequence = protein_sequence
        self.csv_path = csv_path
        self.output_dir = Path(output_dir)

        # Set default thresholds if not provided
        self.target_r2 = target_r2 if target_r2 is not None else 0.6
        self.target_spearman = target_spearman if target_spearman is not None else 0.7
        self.target_roc_auc = target_roc_auc if target_roc_auc is not None else 0.8

        self.max_iterations = max_iterations
        self.use_msa_server = use_msa_server
        self.msa_path = msa_path
        self.roc_threshold = roc_threshold

        # Chain IDs and constraints
        self.protein_chain_id = protein_chain_id
        self.ligand_chain_id = ligand_chain_id
        self.pocket_residues = pocket_residues or []
        self.contact_residues = contact_residues or []
        self.template_force = template_force
        self.template_threshold = template_threshold

        # Create output directories first (before logging setup)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.yaml_dir = self.output_dir / "yaml_inputs"
        self.yaml_dir.mkdir(exist_ok=True)
        self.predictions_dir = self.output_dir / "predictions"
        self.predictions_dir.mkdir(exist_ok=True)

        # Handle template file - convert PDB to CIF if needed
        self.template_file = None
        if template_file:
            template_path = Path(template_file)
            if template_path.suffix.lower() == '.pdb':
                # Convert PDB to CIF
                print(f"Converting PDB template to CIF: {template_file}")
                cif_dir = self.output_dir / "templates"
                cif_dir.mkdir(exist_ok=True)
                self.template_file = convert_pdb_to_cif(template_file, str(cif_dir))
            elif template_path.suffix.lower() == '.cif':
                self.template_file = str(template_path)
            else:
                raise ValueError(f"Template file must be .pdb or .cif: {template_file}")

        # Setup logging (after output_dir is created)
        self._setup_logging(log_level)

        # Load SAR data
        self.sar_data: List[SARData] = []
        self.iteration_results: List[Dict] = []

        self.logger.info(f"Initialized Boltz-2 SAR Iterator")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Convergence thresholds (any one triggers stop):")
        self.logger.info(f"  - R² >= {self.target_r2}")
        self.logger.info(f"  - Spearman R >= {self.target_spearman}")
        self.logger.info(f"  - ROC-AUC >= {self.target_roc_auc}")
        self.logger.info(f"Max iterations: {self.max_iterations}")
        self.logger.info(f"Protein chain: {self.protein_chain_id}, Ligand chain: {self.ligand_chain_id}")
        if self.pocket_residues:
            self.logger.info(f"Pocket constraints: {len(self.pocket_residues)} residues")
        if self.contact_residues:
            self.logger.info(f"Contact constraints: {len(self.contact_residues)} contacts")
        if self.template_file:
            self.logger.info(f"Template file: {self.template_file}")

    def _setup_logging(self, log_level: str):
        """Setup logging configuration."""
        self.logger = logging.getLogger("Boltz2SARIterator")
        self.logger.setLevel(getattr(logging, log_level.upper()))

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(getattr(logging, log_level.upper()))
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        # File handler
        log_file = self.output_dir / "boltz2_sar_iterator.log" if hasattr(self, 'output_dir') else "boltz2_sar_iterator.log"
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)

    def load_sar_data(self) -> None:
        """Load SAR data from CSV file."""
        self.logger.info(f"Loading SAR data from {self.csv_path}")

        try:
            df = pd.read_csv(self.csv_path)
        except Exception as e:
            self.logger.error(f"Failed to load CSV: {e}")
            raise

        # Validate required columns
        required_cols = ['SMILES', 'Activity']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            # Try case-insensitive match
            col_map = {col.lower(): col for col in df.columns}
            if 'smiles' in col_map and 'activity' in col_map:
                df = df.rename(columns={
                    col_map['smiles']: 'SMILES',
                    col_map['activity']: 'Activity'
                })
            else:
                raise ValueError(f"CSV missing required columns: {missing_cols}")

        # Load data
        for idx, row in df.iterrows():
            compound_id = row.get('ID', row.get('id', f"compound_{idx+1}"))
            smiles = row['SMILES']
            activity = float(row['Activity'])

            # Validate SMILES and activity
            if pd.isna(smiles) or smiles == "":
                self.logger.warning(f"Skipping compound {compound_id}: empty SMILES")
                continue
            if pd.isna(activity):
                self.logger.warning(f"Skipping compound {compound_id}: missing activity")
                continue

            self.sar_data.append(SARData(
                compound_id=str(compound_id),
                smiles=str(smiles),
                activity=activity
            ))

        self.logger.info(f"Loaded {len(self.sar_data)} compounds")

        if len(self.sar_data) == 0:
            raise ValueError("No valid compounds found in CSV")

    def create_yaml_input(self, compound: SARData, iteration: int) -> Path:
        """
        Create YAML input file for Boltz-2 prediction.

        Args:
            compound: SARData object
            iteration: Current iteration number

        Returns:
            Path to created YAML file
        """
        yaml_data = {
            'version': 1,
            'sequences': [
                {
                    'protein': {
                        'id': self.protein_chain_id,
                        'sequence': self.protein_sequence
                    }
                },
                {
                    'ligand': {
                        'id': self.ligand_chain_id,
                        'smiles': compound.smiles
                    }
                }
            ]
        }

        # Add MSA - use 'empty' if not provided
        if self.msa_path:
            yaml_data['sequences'][0]['protein']['msa'] = self.msa_path
        else:
            yaml_data['sequences'][0]['protein']['msa'] = 'empty'

        # Add properties for affinity prediction
        yaml_data['properties'] = [
            {
                'affinity': {
                    'binder': self.ligand_chain_id
                }
            }
        ]

        # Add constraints if specified
        constraints = []

        # Add pocket constraints
        if self.pocket_residues:
            pocket_constraint = {
                'pocket': {
                    'binder': self.ligand_chain_id,
                    'contacts': [[chain, residue] for chain, residue in self.pocket_residues]
                }
            }
            constraints.append(pocket_constraint)

        # Add contact constraints
        for contact in self.contact_residues:
            (chain1, res1), (chain2, res2), max_dist = contact
            contact_constraint = {
                'contact': {
                    'token1': [chain1, res1],
                    'token2': [chain2, res2],
                    'max_distance': max_dist,
                    'force': True
                }
            }
            constraints.append(contact_constraint)

        if constraints:
            yaml_data['constraints'] = constraints

        # Add template if specified
        if self.template_file:
            yaml_data['templates'] = [{
                'cif': self.template_file,
                'force': self.template_force,
                'threshold': self.template_threshold
            }]

        yaml_filename = f"iter_{iteration}_{compound.compound_id}.yaml"
        yaml_path = self.yaml_dir / yaml_filename

        with open(yaml_path, 'w') as f:
            yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)

        self.logger.debug(f"Created YAML input: {yaml_path}")
        return yaml_path

    def run_boltz_prediction(self, yaml_path: Path, iteration: int) -> bool:
        """
        Run Boltz-2 prediction for a given YAML input.

        Args:
            yaml_path: Path to YAML input file
            iteration: Current iteration number

        Returns:
            True if prediction succeeded, False otherwise
        """
        cmd = ['boltz', 'predict', str(yaml_path)]

        if self.use_msa_server and not self.msa_path:
            cmd.append('--use_msa_server')

        # Set output directory
        cmd.extend(['--out_dir', str(self.predictions_dir)])

        self.logger.debug(f"Running command: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600  # 10 minute timeout
            )

            if result.returncode == 0:
                self.logger.debug(f"Prediction succeeded for {yaml_path.name}")
                return True
            else:
                self.logger.error(f"Prediction failed for {yaml_path.name}")
                self.logger.error(f"stderr: {result.stderr}")
                return False

        except subprocess.TimeoutExpired:
            self.logger.error(f"Prediction timed out for {yaml_path.name}")
            return False
        except Exception as e:
            self.logger.error(f"Error running prediction: {e}")
            return False

    def extract_affinity(self, compound: SARData, iteration: int) -> Optional[Dict[str, float]]:
        """
        Extract all predicted affinity metrics from Boltz-2 output files.

        Args:
            compound: SARData object
            iteration: Current iteration number

        Returns:
            Dictionary with all affinity metrics or None if extraction failed
        """
        # Boltz-2 creates nested output directory structure:
        # predictions/boltz_results_<yaml_name>/predictions/<yaml_name>/
        yaml_name = f"iter_{iteration}_{compound.compound_id}"

        # Try nested directory structure (Boltz-2 v2.2+)
        nested_prediction_dir = self.predictions_dir / f"boltz_results_{yaml_name}" / "predictions" / yaml_name

        # Fallback to flat directory structure (older versions)
        flat_prediction_dir = self.predictions_dir / yaml_name

        # Determine which directory exists
        if nested_prediction_dir.exists():
            prediction_dir = nested_prediction_dir
        elif flat_prediction_dir.exists():
            prediction_dir = flat_prediction_dir
        else:
            self.logger.warning(f"No prediction directory found for {compound.compound_id}")
            return None

        # Look for affinity JSON files
        affinity_files = list(prediction_dir.glob("affinity_*.json"))

        if not affinity_files:
            self.logger.warning(f"No affinity files found for {compound.compound_id}")
            return None

        # Read the first affinity file (there should typically be one)
        affinity_file = affinity_files[0]
        self.logger.debug(f"Reading affinity from {affinity_file}")

        try:
            with open(affinity_file, 'r') as f:
                affinity_data = json.load(f)

            # Extract all available metrics
            metrics = {}

            # Main predictions
            if 'affinity_pred_value' in affinity_data:
                metrics['affinity_value'] = float(affinity_data['affinity_pred_value'])
            if 'affinity_probability_binary' in affinity_data:
                metrics['probability_binary'] = float(affinity_data['affinity_probability_binary'])

            # Ensemble predictions
            ensemble_values = []
            ensemble_probs = []

            for i in range(1, 10):  # Check for up to 9 ensemble members
                value_key = f'affinity_pred_value{i}'
                prob_key = f'affinity_probability_binary{i}'

                if value_key in affinity_data:
                    ensemble_values.append(float(affinity_data[value_key]))
                if prob_key in affinity_data:
                    ensemble_probs.append(float(affinity_data[prob_key]))

            # Calculate ensemble means if available
            if ensemble_values:
                metrics['ensemble_value_mean'] = np.mean(ensemble_values)
                metrics['ensemble_value_std'] = np.std(ensemble_values)
            if ensemble_probs:
                metrics['ensemble_prob_mean'] = np.mean(ensemble_probs)
                metrics['ensemble_prob_std'] = np.std(ensemble_probs)

            if not metrics:
                self.logger.warning(f"No affinity metrics found in {affinity_file}")
                return None

            self.logger.debug(f"Extracted metrics: {metrics}")
            return metrics

        except Exception as e:
            self.logger.error(f"Error reading affinity file: {e}")
            return None

    def calculate_correlation_metrics(
        self,
        experimental: List[float],
        predicted: List[float],
        roc_threshold: Optional[float] = None
    ) -> Dict[str, float]:
        """
        Calculate multiple correlation metrics between experimental and predicted values.

        Args:
            experimental: List of experimental activity values
            predicted: List of predicted affinity values
            roc_threshold: Optional threshold for ROC-AUC calculation (e.g., IC50 < threshold = binder)

        Returns:
            Dictionary with correlation metrics: r2, spearman, pearson, combined, roc_auc
        """
        exp_array = np.array(experimental)
        pred_array = np.array(predicted)

        # Remove any NaN values
        mask = ~(np.isnan(exp_array) | np.isnan(pred_array))
        exp_array = exp_array[mask]
        pred_array = pred_array[mask]

        if len(exp_array) < 2:
            self.logger.warning("Not enough data points for correlation calculation")
            return {'r2': 0.0, 'spearman': 0.0, 'pearson': 0.0, 'combined': 0.0, 'roc_auc': 0.0}

        # Calculate Pearson correlation (linear)
        pearson_r = np.corrcoef(exp_array, pred_array)[0, 1]
        if np.isnan(pearson_r):
            pearson_r = 0.0

        # Calculate R² (Pearson²) - coefficient of determination for correlation
        r2 = pearson_r ** 2

        # Calculate Spearman correlation (rank-based, better for monotonic relationships)
        spearman_r, _ = spearmanr(exp_array, pred_array)
        if np.isnan(spearman_r):
            spearman_r = 0.0

        # Combined metric: max of R² and Spearman²
        combined = max(r2, spearman_r ** 2)

        # Calculate ROC-AUC if threshold is provided
        roc_auc = 0.0
        if roc_threshold is not None and len(exp_array) >= 2:
            try:
                # Define binary labels based on threshold
                # For IC50: values below threshold are binders (label=1)
                binary_labels = (exp_array < roc_threshold).astype(int)

                # Check if we have both classes
                if len(np.unique(binary_labels)) == 2:
                    # Lower predicted affinity = better binding, so invert for ROC
                    roc_auc = roc_auc_score(binary_labels, -pred_array)
                else:
                    self.logger.debug("ROC-AUC: Only one class present, skipping")
            except Exception as e:
                self.logger.debug(f"ROC-AUC calculation failed: {e}")
                roc_auc = 0.0

        self.logger.debug(
            f"Metrics - R²: {r2:.4f}, Pearson: {pearson_r:.4f}, "
            f"Spearman: {spearman_r:.4f}, Combined: {combined:.4f}"
        )
        if roc_threshold is not None:
            self.logger.debug(f"ROC-AUC (threshold={roc_threshold}): {roc_auc:.4f}")

        return {
            'r2': r2,
            'spearman': spearman_r,
            'pearson': pearson_r,
            'combined': combined,
            'roc_auc': roc_auc
        }

    def calculate_r2(self, experimental: List[float], predicted: List[float]) -> float:
        """
        Legacy method for backward compatibility. Use calculate_correlation_metrics instead.

        Returns the combined metric (max of R² and Spearman²).
        """
        metrics = self.calculate_correlation_metrics(experimental, predicted)
        return metrics['combined']

    def run_iteration(self, iteration: int) -> Tuple[float, int, int, Dict[str, Dict[str, float]]]:
        """
        Run a single iteration of Boltz-2 predictions for all compounds.

        Args:
            iteration: Current iteration number

        Returns:
            Tuple of (combined metric, successful predictions, failed predictions, all metric correlations)
        """
        self.logger.info(f"\n{'='*60}")
        self.logger.info(f"Starting Iteration {iteration}")
        self.logger.info(f"{'='*60}")

        successful = 0
        failed = 0

        for i, compound in enumerate(self.sar_data, 1):
            self.logger.info(f"Processing compound {i}/{len(self.sar_data)}: {compound.compound_id}")

            # Create YAML input
            yaml_path = self.create_yaml_input(compound, iteration)

            # Run prediction
            if self.run_boltz_prediction(yaml_path, iteration):
                # Extract all affinity metrics
                metrics_dict = self.extract_affinity(compound, iteration)
                if metrics_dict is not None:
                    compound.predicted_metrics = metrics_dict
                    # Use affinity_value as primary metric for backward compatibility
                    compound.predicted_affinity = metrics_dict.get('affinity_value', None)
                    successful += 1
                else:
                    failed += 1
            else:
                failed += 1

        # Calculate correlations for all available metrics
        experimental = [c.activity for c in self.sar_data if c.predicted_metrics is not None]

        # Collect all metric types available
        all_metric_names = set()
        for c in self.sar_data:
            if c.predicted_metrics:
                all_metric_names.update(c.predicted_metrics.keys())

        # Calculate correlations for each metric type
        metric_correlations = {}
        best_combined = 0.0
        best_metric_name = None

        self.logger.info(f"\nIteration {iteration} Summary:")
        self.logger.info(f"  Successful predictions: {successful}")
        self.logger.info(f"  Failed predictions: {failed}")

        if len(experimental) > 0:
            self.logger.info(f"\nCorrelations for each metric:")

            for metric_name in sorted(all_metric_names):
                # Get predictions for this metric (only for compounds that have it)
                predicted = []
                exp_subset = []
                for c in self.sar_data:
                    if c.predicted_metrics and metric_name in c.predicted_metrics:
                        predicted.append(c.predicted_metrics[metric_name])
                        exp_subset.append(c.activity)

                if len(predicted) >= 2:
                    corr_metrics = self.calculate_correlation_metrics(exp_subset, predicted, self.roc_threshold)
                    metric_correlations[metric_name] = corr_metrics

                    self.logger.info(f"\n  {metric_name}:")
                    self.logger.info(f"    R² (Pearson²): {corr_metrics['r2']:.4f}")
                    self.logger.info(f"    Spearman R: {corr_metrics['spearman']:.4f}")
                    self.logger.info(f"    Combined: {corr_metrics['combined']:.4f}")
                    if self.roc_threshold is not None and corr_metrics['roc_auc'] > 0:
                        self.logger.info(f"    ROC-AUC: {corr_metrics['roc_auc']:.4f}")

                    # Track best metric
                    if corr_metrics['combined'] > best_combined:
                        best_combined = corr_metrics['combined']
                        best_metric_name = metric_name

            if best_metric_name:
                self.logger.info(f"\n  ✓ Best metric: {best_metric_name} (Combined={best_combined:.4f})")

        # Use best metric's combined score for convergence, or fallback to affinity_value
        if best_metric_name and best_metric_name in metric_correlations:
            metrics = metric_correlations[best_metric_name]
        elif 'affinity_value' in metric_correlations:
            metrics = metric_correlations['affinity_value']
            best_metric_name = 'affinity_value'
        else:
            metrics = {'r2': 0.0, 'spearman': 0.0, 'pearson': 0.0, 'combined': 0.0, 'roc_auc': 0.0}
            best_metric_name = 'none'

        # Store iteration results
        iteration_data = {
            'iteration': iteration,
            'r2': metrics['r2'],
            'spearman': metrics['spearman'],
            'pearson': metrics['pearson'],
            'combined': metrics['combined'],
            'best_metric': best_metric_name,
            'all_metrics': metric_correlations,  # Store correlations for all metrics
            'successful': successful,
            'failed': failed,
            'timestamp': datetime.now().isoformat()
        }
        if self.roc_threshold is not None:
            iteration_data['roc_auc'] = metrics['roc_auc']
            iteration_data['roc_threshold'] = self.roc_threshold

        self.iteration_results.append(iteration_data)

        return metrics['combined'], successful, failed, metric_correlations

    def check_convergence(self, metric_correlations: Dict[str, Dict[str, float]]) -> Tuple[bool, Optional[str], Optional[str]]:
        """
        Check if any metric meets any convergence threshold.

        Args:
            metric_correlations: Dictionary mapping metric names to their correlation results

        Returns:
            Tuple of (converged, metric_name, threshold_name)
        """
        for metric_name, corr_metrics in metric_correlations.items():
            # Check R² threshold
            if corr_metrics['r2'] >= self.target_r2:
                return True, metric_name, f"R² >= {self.target_r2}"

            # Check Spearman threshold
            if abs(corr_metrics['spearman']) >= self.target_spearman:
                return True, metric_name, f"Spearman R >= {self.target_spearman}"

            # Check ROC-AUC threshold (if ROC was calculated)
            if corr_metrics.get('roc_auc', 0) >= self.target_roc_auc:
                return True, metric_name, f"ROC-AUC >= {self.target_roc_auc}"

        return False, None, None

    def run(self) -> Dict:
        """
        Run the iterative Boltz-2 SAR correlation process.

        Returns:
            Dictionary with final results and statistics
        """
        self.logger.info("\n" + "="*60)
        self.logger.info("Starting Boltz-2 SAR Iterator")
        self.logger.info("="*60)

        # Load SAR data
        self.load_sar_data()

        # Run iterations
        converged = False
        final_r2 = 0.0
        convergence_metric = None
        convergence_threshold = None

        for iteration in range(1, self.max_iterations + 1):
            r2, successful, failed, metric_correlations = self.run_iteration(iteration)
            final_r2 = r2

            # Check if any metric meets any threshold
            converged, conv_metric_name, conv_threshold = self.check_convergence(metric_correlations)

            if converged:
                convergence_metric = conv_metric_name
                convergence_threshold = conv_threshold
                self.logger.info(f"\n{'='*60}")
                self.logger.info(f"✓ CONVERGED!")
                self.logger.info(f"  Metric: {conv_metric_name}")
                self.logger.info(f"  Threshold: {conv_threshold}")
                self.logger.info(f"  Value: R²={metric_correlations[conv_metric_name]['r2']:.4f}, "
                               f"Spearman={metric_correlations[conv_metric_name]['spearman']:.4f}")
                if metric_correlations[conv_metric_name].get('roc_auc', 0) > 0:
                    self.logger.info(f"  ROC-AUC: {metric_correlations[conv_metric_name]['roc_auc']:.4f}")
                self.logger.info(f"{'='*60}")
                break

        # Generate final report
        results = self._generate_final_report(converged, final_r2)

        return results

    def _generate_final_report(self, converged: bool, final_r2: float) -> Dict:
        """Generate and save final report."""
        self.logger.info("\n" + "="*60)
        self.logger.info("Final Results")
        self.logger.info("="*60)

        # Create results DataFrame
        results_data = []
        for compound in self.sar_data:
            results_data.append({
                'Compound_ID': compound.compound_id,
                'SMILES': compound.smiles,
                'Experimental_Activity': compound.activity,
                'Predicted_Affinity': compound.predicted_affinity
            })

        df_results = pd.DataFrame(results_data)
        results_csv = self.output_dir / "final_results.csv"
        df_results.to_csv(results_csv, index=False)
        self.logger.info(f"Saved results to {results_csv}")

        # Save iteration history
        df_iterations = pd.DataFrame(self.iteration_results)
        iterations_csv = self.output_dir / "iteration_history.csv"
        df_iterations.to_csv(iterations_csv, index=False)
        self.logger.info(f"Saved iteration history to {iterations_csv}")

        # Summary
        summary = {
            'converged': converged,
            'final_r2': final_r2,
            'target_r2': self.target_r2,
            'iterations_run': len(self.iteration_results),
            'max_iterations': self.max_iterations,
            'total_compounds': len(self.sar_data),
            'successful_predictions': sum(1 for c in self.sar_data if c.predicted_affinity is not None),
            'output_directory': str(self.output_dir)
        }

        # Save summary JSON
        summary_file = self.output_dir / "summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        self.logger.info(f"Saved summary to {summary_file}")

        # Print summary
        self.logger.info(f"\nConverged: {converged}")
        self.logger.info(f"Final R²: {final_r2:.4f}")
        self.logger.info(f"Iterations run: {len(self.iteration_results)}/{self.max_iterations}")
        self.logger.info(f"Successful predictions: {summary['successful_predictions']}/{len(self.sar_data)}")

        return summary


def load_config_file(config_path: str) -> Dict:
    """Load configuration from JSON or YAML file."""
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    if config_path.suffix in ['.yaml', '.yml']:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    elif config_path.suffix == '.json':
        with open(config_path, 'r') as f:
            config = json.load(f)
    else:
        raise ValueError(f"Config file must be .json, .yaml, or .yml: {config_path}")

    return config


def parse_pocket_residues(pocket_str: str) -> List[Tuple[str, int]]:
    """Parse pocket residues from string format: A:107,B:98"""
    residues = []
    for item in pocket_str.split(','):
        chain, res = item.strip().split(':')
        residues.append((chain.strip(), int(res.strip())))
    return residues


def parse_contact_residues(contact_str: str) -> List[Tuple[Tuple[str, int], Tuple[str, int], float]]:
    """Parse contact residues from string format: A:20-B:27:5.0,A:21-B:31:5.0"""
    contacts = []
    for item in contact_str.split(','):
        item = item.strip()
        if not item:
            continue

        # Split by '-' to separate the two residues
        parts = item.split('-')
        if len(parts) != 2:
            continue

        # Parse first residue and chain (e.g., "A:20")
        token1_parts = parts[0].strip().split(':')
        if len(token1_parts) != 2:
            continue
        chain1 = token1_parts[0].strip()
        res1 = int(token1_parts[1].strip())

        # Parse second residue, chain, and distance (e.g., "B:27:5.0")
        token2_parts = parts[1].strip().split(':')
        if len(token2_parts) != 3:
            continue
        chain2 = token2_parts[0].strip()
        res2 = int(token2_parts[1].strip())
        max_dist = float(token2_parts[2].strip())

        contacts.append(((chain1, res1), (chain2, res2), max_dist))

    return contacts


def main():
    """Example usage of Boltz2SARIterator."""
    import argparse

    parser = argparse.ArgumentParser(
        description='Iterative Boltz-2 SAR correlation tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage (command-line):
  python boltz2_sar_iterator.py \\
    --protein MVTPEGNVSLVDESLLVGVTDEDRAV... \\
    --csv compounds.csv \\
    --target-r2 0.7 \\
    --max-iterations 10 \\
    --protein-chain A \\
    --ligand-chain L \\
    --pocket-residues "A:107,A:98" \\
    --template template.cif

Example usage (config file):
  python boltz2_sar_iterator.py --config config.json

CSV file should contain at minimum:
  - SMILES: Chemical structure in SMILES format
  - Activity: Experimental activity values
  - ID (optional): Compound identifier
        """
    )

    # Config file option
    parser.add_argument(
        '--config',
        help='Path to config file (.json or .yaml) with all parameters'
    )

    # Required arguments (if not using config)
    parser.add_argument(
        '--protein',
        help='Protein amino acid sequence'
    )
    parser.add_argument(
        '--csv',
        help='Path to CSV file with SMILES and Activity columns'
    )

    # Optional arguments
    parser.add_argument(
        '--output-dir',
        default='./boltz2_sar_output',
        help='Output directory (default: ./boltz2_sar_output)'
    )
    parser.add_argument(
        '--target-r2',
        type=float,
        default=None,
        help='Target R² threshold - converges if any metric R² >= this (default: 0.6)'
    )
    parser.add_argument(
        '--target-spearman',
        type=float,
        default=None,
        help='Target Spearman R threshold - converges if any metric |Spearman| >= this (default: 0.7)'
    )
    parser.add_argument(
        '--target-roc-auc',
        type=float,
        default=None,
        help='Target ROC-AUC threshold - converges if any metric ROC-AUC >= this (default: 0.8)'
    )
    parser.add_argument(
        '--max-iterations',
        type=int,
        default=10,
        help='Maximum number of iterations (default: 10)'
    )
    parser.add_argument(
        '--msa-path',
        help='Path to pre-computed MSA file (.a3m)'
    )
    parser.add_argument(
        '--no-msa-server',
        action='store_true',
        help='Do not use MSA server'
    )
    parser.add_argument(
        '--protein-chain',
        default='A',
        help='Protein chain ID (default: A)'
    )
    parser.add_argument(
        '--ligand-chain',
        default='L',
        help='Ligand chain ID (default: L)'
    )
    parser.add_argument(
        '--pocket-residues',
        help='Pocket residues as chain:residue,chain:residue (e.g., A:107,A:98)'
    )
    parser.add_argument(
        '--contact-residues',
        help='Contact residues as chain1:res1-chain2:res2:distance,... (e.g., A:20-B:27:5.0)'
    )
    parser.add_argument(
        '--template',
        help='Path to template CIF file'
    )
    parser.add_argument(
        '--template-force',
        action='store_true',
        default=True,
        help='Force template usage (default: True)'
    )
    parser.add_argument(
        '--template-threshold',
        type=int,
        default=2,
        help='Template threshold (default: 2)'
    )
    parser.add_argument(
        '--roc-threshold',
        type=float,
        help='Threshold for ROC-AUC calculation (e.g., 10000 for IC50 < 10μM = binder)'
    )
    parser.add_argument(
        '--log-level',
        default='INFO',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        help='Logging level (default: INFO)'
    )

    args = parser.parse_args()

    # Load configuration
    if args.config:
        print(f"Loading configuration from {args.config}")
        config = load_config_file(args.config)

        # Override with command-line arguments if provided
        protein_sequence = args.protein or config.get('protein_sequence')
        csv_path = args.csv or config.get('csv_path')
        output_dir = args.output_dir if args.output_dir != './boltz2_sar_output' else config.get('output_dir', './boltz2_sar_output')
        target_r2 = args.target_r2 if args.target_r2 != 0.7 else config.get('target_r2', 0.7)
        max_iterations = args.max_iterations if args.max_iterations != 10 else config.get('max_iterations', 10)
        msa_path = args.msa_path or config.get('msa_path')
        use_msa_server = not args.no_msa_server if not args.no_msa_server else config.get('use_msa_server', True)
        protein_chain_id = args.protein_chain if args.protein_chain != 'A' else config.get('protein_chain_id', 'A')
        ligand_chain_id = args.ligand_chain if args.ligand_chain != 'L' else config.get('ligand_chain_id', 'L')

        # Parse pocket and contact residues from config
        pocket_residues = config.get('pocket_residues', [])
        contact_residues = config.get('contact_residues', [])

        template_file = args.template or config.get('template_file')
        template_force = config.get('template_force', True)
        template_threshold = config.get('template_threshold', 2)
        roc_threshold = args.roc_threshold or config.get('roc_threshold')
        log_level = args.log_level if args.log_level != 'INFO' else config.get('log_level', 'INFO')

    else:
        # Use command-line arguments only
        if not args.protein or not args.csv:
            parser.error("--protein and --csv are required when not using --config")

        protein_sequence = args.protein
        csv_path = args.csv
        output_dir = args.output_dir
        target_r2 = args.target_r2
        max_iterations = args.max_iterations
        msa_path = args.msa_path
        use_msa_server = not args.no_msa_server
        protein_chain_id = args.protein_chain
        ligand_chain_id = args.ligand_chain

        # Parse pocket and contact residues from command-line
        pocket_residues = parse_pocket_residues(args.pocket_residues) if args.pocket_residues else []
        contact_residues = parse_contact_residues(args.contact_residues) if args.contact_residues else []

        template_file = args.template
        template_force = args.template_force
        template_threshold = args.template_threshold
        roc_threshold = args.roc_threshold
        log_level = args.log_level

    # Validate arguments
    if not use_msa_server and not msa_path:
        parser.error("--no-msa-server requires --msa-path")

    # Create iterator
    iterator = Boltz2SARIterator(
        protein_sequence=protein_sequence,
        csv_path=csv_path,
        output_dir=output_dir,
        target_r2=config.get('target_r2') if config else args.target_r2,
        target_spearman=config.get('target_spearman') if config else args.target_spearman,
        target_roc_auc=config.get('target_roc_auc') if config else args.target_roc_auc,
        max_iterations=max_iterations,
        use_msa_server=use_msa_server,
        msa_path=msa_path,
        protein_chain_id=protein_chain_id,
        ligand_chain_id=ligand_chain_id,
        pocket_residues=pocket_residues,
        contact_residues=contact_residues,
        template_file=template_file,
        template_force=template_force,
        template_threshold=template_threshold,
        roc_threshold=roc_threshold,
        log_level=log_level
    )

    # Run iterations
    results = iterator.run()

    return results


if __name__ == '__main__':
    main()
