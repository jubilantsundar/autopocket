#!/usr/bin/env python3
"""
Dataset Preparation Script for SAR Workflow Case Study
Prepares curated single-source datasets from BindingDB
Focuses on datasets from single publications to ensure assay consistency
"""

import pandas as pd
import numpy as np
from pathlib import Path
import requests
import gzip
from io import StringIO
from typing import List, Tuple, Dict

# Try to import optional dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not installed. Install with: pip install rdkit")


def download_bindingdb_data(output_file: str = "BindingDB_All.tsv"):
    """
    Download BindingDB full database
    
    Returns:
        Path to downloaded file
    """
    url = "https://www.bindingdb.org/bind/downloads/BindingDB_All_202501_tsv.zip"
    
    print(f"Downloading BindingDB database...")
    print(f"This may take several minutes (file is ~500 MB)...")
    
    # Note: BindingDB updates quarterly, adjust URL if needed
    # Alternative: use wget or curl for large files
    
    print("\nPlease download manually from:")
    print("https://www.bindingdb.org/bind/downloads/")
    print("Select: 'BindingDB_All' TSV format")
    print(f"Save as: {output_file}")
    
    return output_file


def load_bindingdb(file_path: str, target_name: str = None, uniprot_id: str = None) -> pd.DataFrame:
    """
    Load and filter BindingDB data
    
    Args:
        file_path: Path to BindingDB TSV file
        target_name: Target name to filter (e.g., "EGFR")
        uniprot_id: UniProt ID to filter (e.g., "P00533")
    
    Returns:
        Filtered DataFrame
    """
    print(f"Loading BindingDB data from {file_path}...")
    
    # Key columns we need
    usecols = [
        'BindingDB Ligand Name',
        'Ligand SMILES',
        'Target Name',
        'UniProt (SwissProt) Primary ID of Target Chain',
        'Ki (nM)',
        'IC50 (nM)',
        'Kd (nM)',
        'Article DOI',
        'PubMed ID',
        'Institution',
        'Patent Number',
        'BindingDB MonomerID',
    ]
    
    # Load data
    df = pd.read_csv(file_path, sep='\t', usecols=usecols, low_memory=False)
    
    print(f"Loaded {len(df)} total records")
    
    # Filter by target
    if target_name:
        df = df[df['Target Name'].str.contains(target_name, case=False, na=False)]
        print(f"Filtered to {len(df)} records for {target_name}")
    
    if uniprot_id:
        df = df[df['UniProt (SwissProt) Primary ID of Target Chain'] == uniprot_id]
        print(f"Filtered to {len(df)} records for UniProt {uniprot_id}")
    
    return df


def find_single_source_datasets(df: pd.DataFrame, min_compounds: int = 15, 
                                 min_fold_change: float = 3.0) -> List[Dict]:
    """
    Find datasets from single publications with sufficient activity range
    
    Args:
        df: BindingDB DataFrame
        min_compounds: Minimum number of compounds in a dataset
        min_fold_change: Minimum fold change in activity values
    
    Returns:
        List of dataset dictionaries with metadata
    """
    print("\nSearching for single-source datasets...")
    
    datasets = []
    
    # Group by publication (DOI or PubMed ID)
    for source_col in ['Article DOI', 'PubMed ID']:
        grouped = df.groupby(source_col)
        
        for source_id, group in grouped:
            if pd.isna(source_id) or source_id == '':
                continue
            
            # Try each activity type
            for activity_col in ['IC50 (nM)', 'Ki (nM)', 'Kd (nM)']:
                # Filter to rows with this activity
                activity_data = group[pd.notna(group[activity_col])].copy()
                
                if len(activity_data) < min_compounds:
                    continue
                
                # Remove duplicates (same SMILES)
                activity_data = activity_data.drop_duplicates(subset=['Ligand SMILES'])
                
                if len(activity_data) < min_compounds:
                    continue
                
                # Check activity range
                activities = pd.to_numeric(activity_data[activity_col], errors='coerce')
                activities = activities[activities > 0]  # Remove invalid values
                
                if len(activities) < min_compounds:
                    continue
                
                min_act = activities.min()
                max_act = activities.max()
                fold_change = max_act / min_act
                
                if fold_change >= min_fold_change:
                    datasets.append({
                        'source_type': source_col,
                        'source_id': source_id,
                        'activity_type': activity_col.replace(' (nM)', ''),
                        'n_compounds': len(activities),
                        'min_activity': min_act,
                        'max_activity': max_act,
                        'fold_change': fold_change,
                        'data': activity_data,
                        'target_name': activity_data['Target Name'].iloc[0],
                        'uniprot_id': activity_data['UniProt (SwissProt) Primary ID of Target Chain'].iloc[0]
                    })
    
    # Sort by fold change (descending) and number of compounds (descending)
    datasets.sort(key=lambda x: (x['fold_change'], x['n_compounds']), reverse=True)
    
    print(f"\nFound {len(datasets)} single-source datasets meeting criteria")
    
    return datasets


def display_dataset_options(datasets: List[Dict], top_n: int = 10):
    """Display top dataset options"""
    
    print("\n" + "="*100)
    print("TOP SINGLE-SOURCE DATASETS")
    print("="*100)
    print(f"{'#':<4} {'Source':<20} {'Target':<15} {'Activity':<8} {'N':<5} {'Range (nM)':<25} {'Fold':<8}")
    print("-"*100)
    
    for i, ds in enumerate(datasets[:top_n], 1):
        source_str = ds['source_id'][:18] if len(str(ds['source_id'])) > 18 else ds['source_id']
        target_str = ds['target_name'][:13] if len(str(ds['target_name'])) > 13 else ds['target_name']
        range_str = f"{ds['min_activity']:.1f} - {ds['max_activity']:.1f}"
        
        print(f"{i:<4} {source_str:<20} {target_str:<15} {ds['activity_type']:<8} "
              f"{ds['n_compounds']:<5} {range_str:<25} {ds['fold_change']:.1f}x")
    
    print("="*100 + "\n")


def prepare_dataset_from_source(dataset_info: Dict, output_dir: str = "datasets", 
                                 dataset_name: str = None) -> Tuple[pd.DataFrame, str]:
    """
    Prepare a complete dataset from selected source
    
    Args:
        dataset_info: Dataset dictionary from find_single_source_datasets
        output_dir: Output directory
        dataset_name: Custom name for dataset files
    
    Returns:
        Tuple of (DataFrame, dataset_name)
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Generate dataset name if not provided
    if dataset_name is None:
        target_clean = dataset_info['target_name'].replace(' ', '_').replace('/', '_')[:20]
        source_clean = str(dataset_info['source_id']).replace('/', '_').replace('.', '_')[:20]
        dataset_name = f"{target_clean}_{source_clean}"
    
    print("\n" + "="*80)
    print(f"PREPARING DATASET: {dataset_name}")
    print("="*80)
    
    # Extract data
    data = dataset_info['data']
    activity_col = dataset_info['activity_type'] + ' (nM)'
    
    # Create clean dataframe
    df_clean = pd.DataFrame({
        'smiles': data['Ligand SMILES'].values,
        'activity': pd.to_numeric(data[activity_col], errors='coerce').values,
        'compound_id': [f"{dataset_name}_{i+1}" for i in range(len(data))]
    })
    
    # Remove any remaining invalid entries
    df_clean = df_clean.dropna()
    df_clean = df_clean[df_clean['activity'] > 0]
    
    # Sort by activity
    df_clean = df_clean.sort_values('activity').reset_index(drop=True)
    
    # Save ligands CSV
    ligands_file = output_path / f"{dataset_name}_ligands.csv"
    df_clean.to_csv(ligands_file, index=False)
    
    # Print summary
    print(f"\nTarget: {dataset_info['target_name']}")
    print(f"UniProt ID: {dataset_info['uniprot_id']}")
    print(f"Source: {dataset_info['source_type']} = {dataset_info['source_id']}")
    print(f"Activity Type: {dataset_info['activity_type']}")
    print(f"Number of compounds: {len(df_clean)}")
    print(f"Activity range: {df_clean['activity'].min():.2f} - {df_clean['activity'].max():.2f} nM")
    print(f"Fold change: {df_clean['activity'].max() / df_clean['activity'].min():.1f}x")
    print(f"\nSaved ligands to: {ligands_file}")
    
    # Fetch protein sequence
    try:
        uniprot_id = dataset_info['uniprot_id']
        if pd.notna(uniprot_id) and uniprot_id != '':
            sequence = get_protein_sequence(uniprot_id)
            
            # Save protein sequence
            protein_file = output_path / f"{dataset_name}_protein.fasta"
            with open(protein_file, 'w') as f:
                f.write(f">sp|{uniprot_id}|{dataset_info['target_name']}\n")
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + "\n")
            
            print(f"Saved protein sequence to: {protein_file}")
            print(f"Protein length: {len(sequence)} residues")
    except Exception as e:
        print(f"Warning: Could not fetch protein sequence: {e}")
    
    print("="*80 + "\n")
    
    return df_clean, dataset_name


def get_protein_sequence(uniprot_id: str) -> str:
    """Fetch protein sequence from UniProt"""
    print(f"Fetching sequence for {uniprot_id} from UniProt...")
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        lines = response.text.split('\n')
        sequence = ''.join([line for line in lines if not line.startswith('>')])
        return sequence
    else:
        raise ValueError(f"Failed to fetch sequence for {uniprot_id}")


def interactive_dataset_selection(bindingdb_file: str, target_name: str = None, 
                                   uniprot_id: str = None, output_dir: str = "datasets"):
    """
    Interactive workflow to select and prepare datasets
    
    Args:
        bindingdb_file: Path to BindingDB TSV file
        target_name: Optional target name filter
        uniprot_id: Optional UniProt ID filter
        output_dir: Output directory
    """
    # Load BindingDB
    df = load_bindingdb(bindingdb_file, target_name=target_name, uniprot_id=uniprot_id)
    
    # Find single-source datasets
    datasets = find_single_source_datasets(df, min_compounds=15, min_fold_change=3.0)
    
    if len(datasets) == 0:
        print("No suitable single-source datasets found!")
        print("Try relaxing criteria or choosing a different target")
        return
    
    # Display options
    display_dataset_options(datasets, top_n=20)
    
    # Ask user to select
    print("Select datasets to prepare (comma-separated numbers, or 'all' for top 5):")
    print("Example: 1,3,5 or 'all'")
    selection = input("> ").strip()
    
    if selection.lower() == 'all':
        selected_indices = list(range(min(5, len(datasets))))
    else:
        try:
            selected_indices = [int(x.strip()) - 1 for x in selection.split(',')]
        except:
            print("Invalid selection")
            return
    
    # Prepare selected datasets
    prepared_datasets = []
    for idx in selected_indices:
        if 0 <= idx < len(datasets):
            df_prepared, name = prepare_dataset_from_source(
                datasets[idx], 
                output_dir=output_dir
            )
            prepared_datasets.append((name, df_prepared))
    
    # Print final summary
    print("\n" + "="*80)
    print("DATASET PREPARATION COMPLETE")
    print("="*80)
    print(f"\nPrepared {len(prepared_datasets)} dataset(s):\n")
    for name, df in prepared_datasets:
        print(f"  - {name}: {len(df)} compounds")
        print(f"    Run with: python sar_workflow.py {output_dir}/{name}_protein.fasta {output_dir}/{name}_ligands.csv")
    print("\n")


def prepare_recommended_datasets(bindingdb_file: str, output_dir: str = "datasets"):
    """
    Prepare recommended EGFR and GPCR datasets automatically
    """
    print("\n" + "="*80)
    print("PREPARING RECOMMENDED DATASETS")
    print("="*80 + "\n")
    
    # EGFR kinase dataset
    print("Searching for EGFR kinase datasets...")
    df_egfr = load_bindingdb(bindingdb_file, target_name="EGFR", uniprot_id="P00533")
    egfr_datasets = find_single_source_datasets(df_egfr, min_compounds=15, min_fold_change=3.0)
    
    if len(egfr_datasets) > 0:
        print("\nPreparing best EGFR dataset...")
        prepare_dataset_from_source(egfr_datasets[0], output_dir=output_dir, dataset_name="EGFR_kinase")
    else:
        print("No suitable EGFR datasets found")
    
    # A2A receptor dataset
    print("\nSearching for Adenosine A2A receptor datasets...")
    df_a2a = load_bindingdb(bindingdb_file, target_name="Adenosine A2a", uniprot_id="P29274")
    a2a_datasets = find_single_source_datasets(df_a2a, min_compounds=15, min_fold_change=3.0)
    
    if len(a2a_datasets) > 0:
        print("\nPreparing best A2A dataset...")
        prepare_dataset_from_source(a2a_datasets[0], output_dir=output_dir, dataset_name="A2A_receptor")
    else:
        print("No suitable A2A datasets found")


def create_example_dataset(output_dir: str = "datasets"):
    """
    Create example dataset for testing (no external data needed)
    """
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    print("\n" + "="*80)
    print("CREATING EXAMPLE TEST DATASET")
    print("="*80 + "\n")
    
    # Example EGFR-like compounds with realistic IC50 values
    # These are simplified structures for testing
    example_data = {
        'smiles': [
            'Cc1ccc(Nc2ncnc3ccccc23)cc1',  # Simple quinazoline
            'COc1ccc(Nc2ncnc3ccccc23)cc1',  # Methoxy variant
            'Cc1ccc(Nc2ncnc3cc(OC)ccc23)cc1',  # Ring substitution
            'Cc1ccc(Nc2ncnc3cc(OCCOC)ccc23)cc1',  # Larger substituent
            'COc1ccc(Nc2ncnc3cc(OCCOC)ccc23)cc1OC',  # More substituted
            'c1ccc(Nc2ncnc3ccccc23)cc1',  # Unsubstituted
            'Fc1ccc(Nc2ncnc3ccccc23)cc1',  # Fluoro
            'Clc1ccc(Nc2ncnc3ccccc23)cc1',  # Chloro
        ],
        'activity': [5.0, 8.0, 15.0, 35.0, 80.0, 200.0, 450.0, 800.0],  # IC50 in nM, 160-fold range
        'compound_id': [f'example_{i+1}' for i in range(8)]
    }
    
    df = pd.DataFrame(example_data)
    
    # Save ligands
    ligands_file = output_path / "example_ligands.csv"
    df.to_csv(ligands_file, index=False)
    print(f"Saved example ligands to: {ligands_file}")
    
    # Create example protein sequence (EGFR kinase domain fragment)
    # This is a simplified sequence for testing
    example_sequence = (
        "FKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVDNPHVCRLLGIC" 
        "LTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLEDRRLVHRDLAARNVLVKTPQHV"
        "KITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILHRIYTHQSDVWSYGVTVWELMTFGSKPYDGIPAS"
        "EISSILEKGERLPQPPICTIDVYMIMVKCWMIDADSRPKFRELIIEFSKMARDPQRYLVIQGDERMHLPS"
        "PTDSNFYRALMDEEDMDDVVDADEYLIPQQGFFSSPSTSRTPLLSSLSATSNNSTVACIDRNGLQSCPIK"
    )
    
    protein_file = output_path / "example_protein.fasta"
    with open(protein_file, 'w') as f:
        f.write(">Example_Kinase_Domain\n")
        for i in range(0, len(example_sequence), 60):
            f.write(example_sequence[i:i+60] + "\n")
    
    print(f"Saved example protein to: {protein_file}")
    
    print("\n" + "-"*80)
    print("EXAMPLE DATASET SUMMARY")
    print("-"*80)
    print(f"Number of compounds: {len(df)}")
    print(f"Activity range: {df['activity'].min():.2f} - {df['activity'].max():.2f} nM")
    print(f"Fold change: {df['activity'].max() / df['activity'].min():.1f}x")
    print("\nRun with:")
    print(f"  python sar_workflow.py {protein_file} {ligands_file}")
    print("-"*80 + "\n")


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Prepare single-source datasets from BindingDB for SAR workflow"
    )
    parser.add_argument(
        "--bindingdb-file",
        default="BindingDB_All.tsv",
        help="Path to BindingDB TSV file"
    )
    parser.add_argument(
        "--target",
        help="Target name to filter (e.g., 'EGFR', 'Adenosine A2a')"
    )
    parser.add_argument(
        "--uniprot",
        help="UniProt ID to filter (e.g., 'P00533' for EGFR)"
    )
    parser.add_argument(
        "--mode",
        choices=["interactive", "recommended", "example"],
        default="interactive",
        help="Mode: interactive selection, recommended datasets, or example dataset"
    )
    parser.add_argument(
        "--output-dir",
        default="datasets",
        help="Output directory (default: datasets)"
    )
    
    args = parser.parse_args()
    
    if args.mode == "example":
        create_example_dataset(output_dir=args.output_dir)
    
    elif args.mode == "recommended":
        if not Path(args.bindingdb_file).exists():
            print(f"Error: BindingDB file not found: {args.bindingdb_file}")
            print("\nPlease download from: https://www.bindingdb.org/bind/downloads/")
            print("Select: 'BindingDB_All' in TSV format")
            return
        
        prepare_recommended_datasets(args.bindingdb_file, output_dir=args.output_dir)
    
    else:  # interactive
        if not Path(args.bindingdb_file).exists():
            print(f"Error: BindingDB file not found: {args.bindingdb_file}")
            print("\nPlease download from: https://www.bindingdb.org/bind/downloads/")
            print("Select: 'BindingDB_All' in TSV format")
            print("\nOr use --mode example to create a test dataset")
            return
        
        interactive_dataset_selection(
            args.bindingdb_file,
            target_name=args.target,
            uniprot_id=args.uniprot,
            output_dir=args.output_dir
        )


if __name__ == "__main__":
    main()
