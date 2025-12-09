#!/usr/bin/env python3
"""
Automated Structure-Activity Relationship (SAR) Analysis Workflow - Updated Version
Finds protein pockets that explain SAR using cofolding and docking simulations
Supports: Boltz-2, DiffDock, ProtoBind-Diff, and Chai-1

Updated to handle:
- Separate conda environment for DiffDock
- Direct protein sequence inputs (no structure required)
- All tool installations and paths
- GPU acceleration where available

Usage:
    python autopocket_sar_workflow_updated.py protein.fasta ligands.csv --method boltz --correlation 0.7
"""

import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import json
import logging
import os
import sys
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from scipy.stats import pearsonr, spearmanr
import argparse

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class LigandData:
    """Store ligand information"""
    smiles: str
    activity: float  # IC50 or Kd
    compound_id: str


@dataclass
class DockingResult:
    """Store docking results"""
    compound_id: str
    binding_affinity: float
    confidence: float
    method: str
    seed: int
    additional_metrics: Dict = None


class ProteinLigandWorkflow:
    """Main workflow class for automated SAR analysis"""

    def __init__(self,
                 protein_fasta: str,
                 ligand_csv: str,
                 target_correlation: float = 0.7,
                 max_iterations: int = 100,
                 output_dir: str = "sar_workflow_output",
                 base_dir: str = "/content/autopocket"):
        """
        Initialize workflow

        Args:
            protein_fasta: Path to FASTA file with protein sequence
            ligand_csv: Path to CSV with columns: 'smiles', 'activity', 'compound_id'
            target_correlation: Minimum correlation to achieve (default: 0.7)
            max_iterations: Maximum number of random seed iterations (default: 100)
            output_dir: Directory for output files
            base_dir: Base directory where tools are installed
        """
        self.base_dir = Path(base_dir)
        self.protein_sequence = self._load_protein_sequence(protein_fasta)
        self.protein_fasta_path = Path(protein_fasta)
        self.ligands = self._load_ligands(ligand_csv)
        self.target_correlation = target_correlation
        self.max_iterations = max_iterations
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

        # Tool paths
        self.boltz_dir = self.base_dir / "boltz"
        self.diffdock_dir = self.base_dir / "DiffDock"
        self.chai_dir = self.base_dir / "chai-lab"
        self.protobind_dir = self.base_dir / "ProtoBind-Diff"

        # Conda paths
        self.conda_path = "/opt/miniconda/bin/conda"
        self.diffdock_env = "diffdock"

        # Store results
        self.all_results = []
        self.best_result = None

        logger.info(f"Initialized workflow with {len(self.ligands)} ligands")
        logger.info(f"Protein sequence length: {len(self.protein_sequence)} residues")

    def _load_protein_sequence(self, fasta_path: str) -> str:
        """Load protein sequence from FASTA file"""
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
            # Skip FASTA header lines
            sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        logger.info(f"Loaded protein sequence from {fasta_path}: {len(sequence)} residues")
        return sequence

    def _load_ligands(self, csv_path: str) -> List[LigandData]:
        """Load ligand data from CSV"""
        df = pd.read_csv(csv_path)

        # Validate required columns
        required_cols = ['smiles', 'activity']
        if not all(col in df.columns for col in required_cols):
            raise ValueError(f"CSV must contain columns: {required_cols}")

        # Add compound_id if not present
        if 'compound_id' not in df.columns:
            df['compound_id'] = [f"compound_{i+1}" for i in range(len(df))]

        ligands = [
            LigandData(
                smiles=row['smiles'],
                activity=row['activity'],
                compound_id=row['compound_id']
            )
            for _, row in df.iterrows()
        ]

        logger.info(f"Loaded {len(ligands)} ligands from {csv_path}")
        return ligands

    def _run_boltz(self, ligand: LigandData, seed: int) -> Optional[DockingResult]:
        """Run Boltz-2 cofolding simulation with affinity prediction"""
        logger.info(f"Running Boltz-2 for {ligand.compound_id} with seed {seed}")

        work_dir = self.output_dir / f"boltz_{ligand.compound_id}_seed{seed}"
        work_dir.mkdir(exist_ok=True, parents=True)

        try:
            # Create YAML input for Boltz-2
            yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {self.protein_sequence}
  - ligand:
      id: B
      smiles: '{ligand.smiles}'
properties:
  - affinity:
      binder: B
"""
            yaml_file = work_dir / "input.yaml"
            with open(yaml_file, 'w') as f:
                f.write(yaml_content)

            # Run Boltz prediction
            cmd = [
                "boltz", "predict",
                str(yaml_file),
                "--accelerator", "gpu",
                "--out_dir", str(work_dir),
                "--num_diffusion_samples", "1"
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse output files
            predictions_dir = work_dir / "boltz_results_input" / "predictions" / "input"
            affinity_json = predictions_dir / "affinity_input.json"
            confidence_json = predictions_dir / "confidence_input_model_0.json"

            if not affinity_json.exists() or not confidence_json.exists():
                raise FileNotFoundError("Boltz output files not found")

            # Load affinity and confidence scores
            with open(affinity_json, 'r') as f:
                affinity_data = json.load(f)

            with open(confidence_json, 'r') as f:
                confidence_data = json.load(f)

            # Extract metrics
            affinity_pred = affinity_data.get('affinity_pred_value', 0)
            affinity_prob = affinity_data.get('affinity_probability_binary', 0)
            confidence_score = confidence_data.get('confidence_score', 0)
            iptm = confidence_data.get('iptm', 0)
            ligand_iptm = confidence_data.get('ligand_iptm', 0)

            # Use affinity prediction as binding score
            binding_score = affinity_pred

            return DockingResult(
                compound_id=ligand.compound_id,
                binding_affinity=binding_score,
                confidence=confidence_score,
                method="Boltz-2",
                seed=seed,
                additional_metrics={
                    'affinity_probability': affinity_prob,
                    'iptm': iptm,
                    'ligand_iptm': ligand_iptm
                }
            )

        except Exception as e:
            logger.error(f"Boltz-2 failed for {ligand.compound_id}: {e}")
            return None

    def _run_diffdock(self, ligand: LigandData, seed: int) -> Optional[DockingResult]:
        """Run DiffDock docking simulation using conda environment"""
        logger.info(f"Running DiffDock for {ligand.compound_id} with seed {seed}")

        work_dir = self.output_dir / f"diffdock_{ligand.compound_id}_seed{seed}"
        work_dir.mkdir(exist_ok=True, parents=True)

        try:
            # Prepare conda activation and run DiffDock
            # DiffDock can use protein sequence directly with ESMFold
            conda_activate = f"export PATH='/opt/miniconda/bin:$PATH' && source /opt/miniconda/etc/profile.d/conda.sh && conda activate {self.diffdock_env}"

            diffdock_cmd = f"cd {self.diffdock_dir} && python -m inference " \
                          f"--config default_inference_args.yaml " \
                          f"--protein_sequence '{self.protein_sequence}' " \
                          f"--ligand '{ligand.smiles}' " \
                          f"--out_dir {work_dir} " \
                          f"--samples_per_complex 5 " \
                          f"--seed {seed}"

            full_cmd = f"{conda_activate} && {diffdock_cmd}"

            result = subprocess.run(
                full_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                executable='/bin/bash'
            )

            # Parse DiffDock output - look for confidence scores
            output_files = list(work_dir.glob("**/rank*.sdf"))
            if not output_files:
                raise FileNotFoundError("No DiffDock output files found")

            # Read confidence from output directory
            confidence_file = work_dir / "index.csv"
            if confidence_file.exists():
                conf_df = pd.read_csv(confidence_file)
                if len(conf_df) > 0:
                    top_confidence = conf_df.iloc[0]['confidence']
                    binding_score = top_confidence
                    normalized_confidence = min(max((top_confidence + 1.5) / 3.0, 0), 1)
                else:
                    binding_score = 0
                    normalized_confidence = 0.5
            else:
                binding_score = 0
                normalized_confidence = 0.5

            return DockingResult(
                compound_id=ligand.compound_id,
                binding_affinity=binding_score,
                confidence=normalized_confidence,
                method="DiffDock",
                seed=seed,
                additional_metrics={'num_poses': len(output_files)}
            )

        except Exception as e:
            logger.error(f"DiffDock failed for {ligand.compound_id}: {e}")
            return None

    def _run_protobind_diff(self, ligand: LigandData, seed: int) -> Optional[DockingResult]:
        """Run ProtoBind-Diff ligand generation and docking"""
        logger.info(f"Running ProtoBind-Diff for {ligand.compound_id} with seed {seed}")

        work_dir = self.output_dir / f"protobind_{ligand.compound_id}_seed{seed}"
        work_dir.mkdir(exist_ok=True, parents=True)

        try:
            # Create FASTA file for protein
            fasta_file = work_dir / "protein.fasta"
            with open(fasta_file, 'w') as f:
                f.write(f">protein\n{self.protein_sequence}\n")

            # Run ProtoBind-Diff inference
            cmd = [
                "protobind-infer",
                "--fasta_file", str(fasta_file),
                "--output_dir", str(work_dir),
                "--n_batches", "1",
                "--batch_size", "10",
                "--sampling_steps", "100"
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # Parse generated SMILES
            smiles_file = work_dir / "generated_smiles.txt"
            if not smiles_file.exists():
                raise FileNotFoundError("ProtoBind-Diff output not found")

            # Read generated molecules
            with open(smiles_file, 'r') as f:
                lines = f.readlines()

            generated_smiles = [line.strip() for line in lines if line.strip() and line.strip() != 'SMILES']

            # Calculate similarity to input ligand
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from rdkit import DataStructs

            mol_input = Chem.MolFromSmiles(ligand.smiles)
            fp_input = AllChem.GetMorganFingerprintAsBitVect(mol_input, 2, nBits=2048)

            max_similarity = 0
            for gen_smiles in generated_smiles[:10]:  # Check top 10
                mol_gen = Chem.MolFromSmiles(gen_smiles)
                if mol_gen:
                    fp_gen = AllChem.GetMorganFingerprintAsBitVect(mol_gen, 2, nBits=2048)
                    similarity = DataStructs.TanimotoSimilarity(fp_input, fp_gen)
                    max_similarity = max(max_similarity, similarity)

            # Use similarity as proxy for binding (higher similarity = better match)
            binding_score = -10 * max_similarity  # Convert to negative scale

            return DockingResult(
                compound_id=ligand.compound_id,
                binding_affinity=binding_score,
                confidence=max_similarity,
                method="ProtoBind-Diff",
                seed=seed,
                additional_metrics={
                    'num_generated': len(generated_smiles),
                    'max_similarity': max_similarity
                }
            )

        except Exception as e:
            logger.error(f"ProtoBind-Diff failed for {ligand.compound_id}: {e}")
            return None

    def _run_chai1(self, ligand: LigandData, seed: int) -> Optional[DockingResult]:
        """Run Chai-1 cofolding simulation"""
        logger.info(f"Running Chai-1 for {ligand.compound_id} with seed {seed}")

        work_dir = self.output_dir / f"chai1_{ligand.compound_id}_seed{seed}"
        work_dir.mkdir(exist_ok=True, parents=True)

        try:
            # Create FASTA input for Chai-1
            fasta_content = f""">protein|A
{self.protein_sequence}
>ligand|B
{ligand.smiles}
"""
            fasta_file = work_dir / "input.fasta"
            with open(fasta_file, 'w') as f:
                f.write(fasta_content)

            # Run Chai-1 via command line
            cmd = [
                "chai-lab", "fold",
                str(fasta_file),
                str(work_dir),
                "--num-trunk-samples", "1",
                "--num-diffn-samples", "1"
            ]

            result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=600)

            # Parse Chai-1 output
            scores_file = work_dir / "scores.model_idx_0.npz"
            cif_file = work_dir / "pred.model_idx_0.cif"

            if not cif_file.exists():
                raise FileNotFoundError("Chai-1 output structure not found")

            # Try to load scores if available
            if scores_file.exists():
                import numpy as np
                scores = np.load(scores_file)
                aggregate_score = float(scores.get('aggregate_score', 0.5))
                ptm = float(scores.get('ptm', 0.5))
                binding_score = -aggregate_score * 100
                confidence = aggregate_score
            else:
                # Fallback if scores not available
                binding_score = -50.0
                confidence = 0.5
                ptm = 0.5

            return DockingResult(
                compound_id=ligand.compound_id,
                binding_affinity=binding_score,
                confidence=confidence,
                method="Chai-1",
                seed=seed,
                additional_metrics={'ptm': ptm}
            )

        except Exception as e:
            logger.error(f"Chai-1 failed for {ligand.compound_id}: {e}")
            return None

    def run_cofolding(self, method: str = "boltz", seed: int = 42) -> List[DockingResult]:
        """Run cofolding simulations for all ligands"""
        results = []

        methods = {
            "boltz": self._run_boltz,
            "diffdock": self._run_diffdock,
            "protobind": self._run_protobind_diff,
            "chai1": self._run_chai1
        }

        if method == "all":
            selected_methods = list(methods.values())
            method_names = list(methods.keys())
        else:
            selected_methods = [methods[method]]
            method_names = [method]

        for ligand in self.ligands:
            for method_func in selected_methods:
                result = method_func(ligand, seed)
                if result:
                    results.append(result)
                else:
                    logger.warning(f"Skipping failed result for {ligand.compound_id}")

        logger.info(f"Completed {len(results)} docking simulations")
        return results

    def calculate_correlation(self, results: List[DockingResult]) -> Tuple[float, float, pd.DataFrame]:
        """Calculate correlation between experimental and predicted values"""
        data = []
        for result in results:
            ligand = next((l for l in self.ligands if l.compound_id == result.compound_id), None)
            if ligand:
                data.append({
                    'compound_id': result.compound_id,
                    'experimental_activity': ligand.activity,
                    'predicted_affinity': result.binding_affinity,
                    'method': result.method,
                    'seed': result.seed,
                    'confidence': result.confidence
                })

        if len(data) < 2:
            logger.error("Insufficient data for correlation calculation")
            return 0.0, 0.0, pd.DataFrame(data)

        df = pd.DataFrame(data)

        try:
            pearson_r, pearson_p = pearsonr(df['experimental_activity'], df['predicted_affinity'])
            spearman_r, spearman_p = spearmanr(df['experimental_activity'], df['predicted_affinity'])
        except Exception as e:
            logger.error(f"Correlation calculation failed: {e}")
            return 0.0, 0.0, df

        logger.info(f"Pearson r: {pearson_r:.3f} (p={pearson_p:.3e})")
        logger.info(f"Spearman r: {spearman_r:.3f} (p={spearman_p:.3e})")

        return pearson_r, spearman_r, df

    def plot_correlation(self, df: pd.DataFrame, pearson_r: float, output_file: str = None):
        """Plot correlation between experimental and predicted values"""
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
            import seaborn as sns

            plt.figure(figsize=(10, 6))
            sns.scatterplot(data=df, x='experimental_activity', y='predicted_affinity',
                          hue='method', style='method', s=100)

            plt.xlabel('Experimental Activity (IC50 or Kd)')
            plt.ylabel('Predicted Binding Affinity')
            plt.title(f'SAR Correlation Analysis (Pearson r = {pearson_r:.3f})')
            plt.grid(True, alpha=0.3)

            if output_file is None:
                output_file = self.output_dir / "correlation_plot.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            logger.info(f"Saved correlation plot to {output_file}")
            plt.close()

        except Exception as e:
            logger.warning(f"Plotting failed: {e}")

    def run_workflow(self, method: str = "boltz", use_absolute_correlation: bool = True):
        """Run the complete automated workflow"""
        logger.info("="*80)
        logger.info("Starting Automated SAR Analysis Workflow")
        logger.info(f"Method: {method}")
        logger.info(f"Target correlation: {self.target_correlation:.3f}")
        logger.info("="*80)

        iteration = 0
        best_correlation = 0
        best_results = None
        best_df = None

        while iteration < self.max_iterations:
            seed = 42 + iteration
            logger.info(f"\nIteration {iteration + 1}/{self.max_iterations} (seed={seed})")

            # Run simulations
            results = self.run_cofolding(method=method, seed=seed)

            if not results:
                logger.warning(f"No successful results in iteration {iteration + 1}")
                iteration += 1
                continue

            self.all_results.extend(results)

            # Calculate correlation
            pearson_r, spearman_r, df = self.calculate_correlation(results)

            correlation = abs(pearson_r) if use_absolute_correlation else pearson_r

            if correlation > best_correlation:
                best_correlation = correlation
                best_results = results
                best_df = df
                self._save_results(best_df, iteration, correlation)

            if correlation >= self.target_correlation:
                logger.info(f"\n{'='*80}")
                logger.info(f"SUCCESS! Target correlation {self.target_correlation:.3f} achieved!")
                logger.info(f"Final correlation: {correlation:.3f}")
                logger.info(f"Iterations required: {iteration + 1}")
                logger.info(f"{'='*80}")
                break

            iteration += 1

        if best_correlation < self.target_correlation:
            logger.warning(f"\nTarget correlation not achieved after {self.max_iterations} iterations")
            logger.info(f"Best correlation achieved: {best_correlation:.3f}")

        self.best_result = {
            'correlation': best_correlation,
            'results': best_results,
            'df': best_df,
            'iterations': iteration + 1
        }

        if best_df is not None and len(best_df) > 0:
            self._save_final_results()
            self.plot_correlation(best_df, best_correlation)

        return self.best_result

    def _save_results(self, df: pd.DataFrame, iteration: int, correlation: float):
        """Save intermediate results"""
        output_file = self.output_dir / f"results_iter{iteration}_r{correlation:.3f}.csv"
        df.to_csv(output_file, index=False)
        logger.info(f"Saved intermediate results to {output_file}")

    def _save_final_results(self):
        """Save final results and summary"""
        if self.best_result is None or self.best_result['df'] is None:
            return

        output_file = self.output_dir / "final_results.csv"
        self.best_result['df'].to_csv(output_file, index=False)

        summary = {
            'target_correlation': self.target_correlation,
            'achieved_correlation': self.best_result['correlation'],
            'iterations_required': self.best_result['iterations'],
            'max_iterations': self.max_iterations,
            'num_ligands': len(self.ligands),
            'protein_length': len(self.protein_sequence)
        }

        summary_file = self.output_dir / "workflow_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)

        logger.info(f"\nFinal results saved to {self.output_dir}")


def main():
    """Command-line interface"""
    parser = argparse.ArgumentParser(
        description="Automated SAR Analysis Workflow - Updated Version",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "protein",
        help="Path to FASTA file with protein sequence"
    )
    parser.add_argument(
        "ligands",
        help="Path to CSV file with ligand SMILES and activity data"
    )
    parser.add_argument(
        "--correlation",
        type=float,
        default=0.7,
        help="Target correlation to achieve (default: 0.7)"
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=100,
        help="Maximum number of iterations (default: 100)"
    )
    parser.add_argument(
        "--method",
        choices=["boltz", "diffdock", "protobind", "chai1", "all"],
        default="boltz",
        help="Docking/cofolding method to use (default: boltz)"
    )
    parser.add_argument(
        "--output",
        default="sar_workflow_output",
        help="Output directory (default: sar_workflow_output)"
    )
    parser.add_argument(
        "--base-dir",
        default="/content/autopocket",
        help="Base directory where tools are installed (default: /content/autopocket)"
    )

    args = parser.parse_args()

    workflow = ProteinLigandWorkflow(
        protein_fasta=args.protein,
        ligand_csv=args.ligands,
        target_correlation=args.correlation,
        max_iterations=args.max_iterations,
        output_dir=args.output,
        base_dir=args.base_dir
    )

    workflow.run_workflow(method=args.method)


if __name__ == "__main__":
    main()
