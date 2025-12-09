#!/usr/bin/env python3
"""
Simple script to extract ligands from any TSV file
User specifies which columns to use for SMILES and activity
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse


def list_columns(tsv_file: str):
    """
    List all columns in the TSV file
    """
    print(f"\nReading columns from: {tsv_file}")
    print("="*80)
    
    df = pd.read_csv(tsv_file, sep='\t', nrows=5)
    
    print(f"\nFound {len(df.columns)} columns:\n")
    
    for i, col in enumerate(df.columns, 1):
        # Show first non-null value as example
        example = df[col].dropna().iloc[0] if len(df[col].dropna()) > 0 else "N/A"
        print(f"{i:3d}. {col}")
        print(f"     Example: {str(example)[:60]}")
        print()
    
    print("="*80)


def extract_ligands(tsv_file: str, smiles_col: str, activity_col: str, 
                   compound_id_col: str = None, output_file: str = "ligands.csv",
                   min_activity: float = None, max_activity: float = None,
                   remove_duplicates: bool = True):
    """
    Extract ligands from TSV file
    
    Args:
        tsv_file: Path to input TSV file
        smiles_col: Column name for SMILES
        activity_col: Column name for activity (IC50, Ki, Kd, etc.)
        compound_id_col: Column name for compound ID (optional)
        output_file: Output CSV file name
        min_activity: Minimum activity value to include (optional)
        max_activity: Maximum activity value to include (optional)
        remove_duplicates: Remove duplicate SMILES (default: True)
    """
    print(f"\nReading data from: {tsv_file}")
    print("="*80)
    
    # Read TSV file
    df = pd.read_csv(tsv_file, sep='\t')
    print(f"Loaded {len(df)} total rows")
    
    # Check if columns exist
    if smiles_col not in df.columns:
        print(f"Error: Column '{smiles_col}' not found!")
        print("Available columns:")
        for col in df.columns:
            print(f"  - {col}")
        return
    
    if activity_col not in df.columns:
        print(f"Error: Column '{activity_col}' not found!")
        print("Available columns:")
        for col in df.columns:
            print(f"  - {col}")
        return
    
    # Extract relevant columns
    if compound_id_col and compound_id_col in df.columns:
        df_extract = df[[smiles_col, activity_col, compound_id_col]].copy()
    else:
        df_extract = df[[smiles_col, activity_col]].copy()
    
    # Rename columns
    df_extract.columns = ['smiles', 'activity'] + ((['compound_id'] if compound_id_col else []))
    
    # Remove rows with missing data
    initial_count = len(df_extract)
    df_extract = df_extract.dropna(subset=['smiles', 'activity'])
    print(f"Removed {initial_count - len(df_extract)} rows with missing SMILES or activity")
    
    # Convert activity to numeric
    df_extract['activity'] = pd.to_numeric(df_extract['activity'], errors='coerce')
    df_extract = df_extract.dropna(subset=['activity'])
    
    # Remove invalid activities (negative or zero)
    df_extract = df_extract[df_extract['activity'] > 0]
    print(f"Kept {len(df_extract)} rows with valid numeric activity values")
    
    # Filter by activity range if specified
    if min_activity is not None:
        df_extract = df_extract[df_extract['activity'] >= min_activity]
        print(f"Filtered to activity >= {min_activity}: {len(df_extract)} rows")
    
    if max_activity is not None:
        df_extract = df_extract[df_extract['activity'] <= max_activity]
        print(f"Filtered to activity <= {max_activity}: {len(df_extract)} rows")
    
    # Remove duplicates based on SMILES
    if remove_duplicates:
        before_dedup = len(df_extract)
        df_extract = df_extract.drop_duplicates(subset=['smiles'], keep='first')
        print(f"Removed {before_dedup - len(df_extract)} duplicate SMILES")
    
    # Add compound IDs if not provided
    if 'compound_id' not in df_extract.columns:
        df_extract['compound_id'] = [f"compound_{i+1}" for i in range(len(df_extract))]
    
    # Sort by activity (lowest to highest)
    df_extract = df_extract.sort_values('activity').reset_index(drop=True)
    
    # Calculate statistics
    min_act = df_extract['activity'].min()
    max_act = df_extract['activity'].max()
    fold_change = max_act / min_act
    
    print("\n" + "="*80)
    print("DATASET SUMMARY")
    print("="*80)
    print(f"Number of compounds: {len(df_extract)}")
    print(f"Activity range: {min_act:.2f} - {max_act:.2f}")
    print(f"Fold change: {fold_change:.1f}x")
    print(f"Mean activity: {df_extract['activity'].mean():.2f}")
    print(f"Median activity: {df_extract['activity'].median():.2f}")
    
    # Check for good SAR range
    if fold_change < 3:
        print(f"\n⚠️  WARNING: Fold change ({fold_change:.1f}x) is less than 3x")
        print("   This may not provide sufficient activity range for SAR analysis")
    else:
        print(f"\n✓ Good activity range for SAR analysis ({fold_change:.1f}x fold change)")
    
    # Save to CSV
    output_path = Path(output_file)
    df_extract[['smiles', 'activity', 'compound_id']].to_csv(output_path, index=False)
    
    print(f"\nSaved {len(df_extract)} compounds to: {output_file}")
    print("="*80 + "\n")
    
    # Show first few entries
    print("First 5 compounds:")
    print(df_extract[['compound_id', 'activity', 'smiles']].head())
    print()
    
    return df_extract


def interactive_mode(tsv_file: str, output_file: str = "ligands.csv"):
    """
    Interactive mode - ask user for column names
    """
    print("\n" + "="*80)
    print("INTERACTIVE LIGAND EXTRACTION")
    print("="*80)
    
    # Show available columns
    list_columns(tsv_file)
    
    # Ask for column names
    print("\nEnter column names (or column numbers):")
    print("-"*80)
    
    smiles_input = input("SMILES column (name or number): ").strip()
    activity_input = input("Activity column (name or number): ").strip()
    compound_id_input = input("Compound ID column (optional, press Enter to skip): ").strip()
    
    # Load dataframe to get column names
    df = pd.read_csv(tsv_file, sep='\t', nrows=1)
    columns = list(df.columns)
    
    # Convert to column names if numbers provided
    def get_column_name(input_str, columns):
        if input_str.isdigit():
            idx = int(input_str) - 1
            if 0 <= idx < len(columns):
                return columns[idx]
            else:
                print(f"Error: Column number {input_str} out of range")
                return None
        else:
            return input_str
    
    smiles_col = get_column_name(smiles_input, columns)
    activity_col = get_column_name(activity_input, columns)
    compound_id_col = get_column_name(compound_id_input, columns) if compound_id_input else None
    
    if not smiles_col or not activity_col:
        print("Error: Invalid column selection")
        return
    
    # Ask for filters
    print("\nOptional filters:")
    print("-"*80)
    min_activity = input("Minimum activity value (press Enter to skip): ").strip()
    max_activity = input("Maximum activity value (press Enter to skip): ").strip()
    
    min_act = float(min_activity) if min_activity else None
    max_act = float(max_activity) if max_activity else None
    
    # Extract ligands
    extract_ligands(
        tsv_file=tsv_file,
        smiles_col=smiles_col,
        activity_col=activity_col,
        compound_id_col=compound_id_col,
        output_file=output_file,
        min_activity=min_act,
        max_activity=max_act,
        remove_duplicates=True
    )


def main():
    parser = argparse.ArgumentParser(
        description="Extract ligands from TSV file for SAR workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Interactive mode - script will ask for column names
  python extract_ligands.py data.tsv --interactive

  # List available columns
  python extract_ligands.py data.tsv --list-columns

  # Direct extraction with column names
  python extract_ligands.py data.tsv \\
      --smiles-col "Ligand SMILES" \\
      --activity-col "IC50 (nM)" \\
      --output ligands.csv

  # With activity range filter
  python extract_ligands.py data.tsv \\
      --smiles-col "SMILES" \\
      --activity-col "Ki (nM)" \\
      --min-activity 1 \\
      --max-activity 1000 \\
      --output egfr_ligands.csv

  # With compound ID column
  python extract_ligands.py data.tsv \\
      --smiles-col "SMILES" \\
      --activity-col "IC50 (nM)" \\
      --compound-id-col "Compound Name" \\
      --output ligands.csv
        """
    )
    
    parser.add_argument(
        "tsv_file",
        help="Path to input TSV file"
    )
    parser.add_argument(
        "--list-columns",
        action="store_true",
        help="List all columns in the TSV file and exit"
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Interactive mode - ask for column names"
    )
    parser.add_argument(
        "--smiles-col",
        help="Name of column containing SMILES strings"
    )
    parser.add_argument(
        "--activity-col",
        help="Name of column containing activity values (IC50, Ki, Kd, etc.)"
    )
    parser.add_argument(
        "--compound-id-col",
        help="Name of column containing compound IDs (optional)"
    )
    parser.add_argument(
        "--output",
        default="ligands.csv",
        help="Output CSV file (default: ligands.csv)"
    )
    parser.add_argument(
        "--min-activity",
        type=float,
        help="Minimum activity value to include"
    )
    parser.add_argument(
        "--max-activity",
        type=float,
        help="Maximum activity value to include"
    )
    parser.add_argument(
        "--keep-duplicates",
        action="store_true",
        help="Keep duplicate SMILES (default: remove duplicates)"
    )
    
    args = parser.parse_args()
    
    # Check if file exists
    if not Path(args.tsv_file).exists():
        print(f"Error: File not found: {args.tsv_file}")
        return
    
    # List columns mode
    if args.list_columns:
        list_columns(args.tsv_file)
        return
    
    # Interactive mode
    if args.interactive:
        interactive_mode(args.tsv_file, output_file=args.output)
        return
    
    # Direct mode - need column names
    if not args.smiles_col or not args.activity_col:
        print("Error: --smiles-col and --activity-col are required")
        print("Use --interactive mode or --list-columns to see available columns")
        return
    
    # Extract ligands
    extract_ligands(
        tsv_file=args.tsv_file,
        smiles_col=args.smiles_col,
        activity_col=args.activity_col,
        compound_id_col=args.compound_id_col,
        output_file=args.output,
        min_activity=args.min_activity,
        max_activity=args.max_activity,
        remove_duplicates=not args.keep_duplicates
    )


if __name__ == "__main__":
    main()
