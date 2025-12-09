#!/usr/bin/env python3
"""
Example usage script for Boltz-2 SAR Iterator.

This demonstrates how to use the tool both as a command-line tool
and programmatically via the Python API.
"""

from boltz2_sar_iterator import Boltz2SARIterator


def example_basic_usage():
    """
    Basic example with minimal parameters.
    """
    print("="*60)
    print("Example 1: Basic Usage")
    print("="*60)

    # Example protein sequence (truncated for demo)
    protein_sequence = (
        "STNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIK"
        "KRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEETEIMIVQAKGRGRGRK"
    )

    iterator = Boltz2SARIterator(
        protein_sequence=protein_sequence,
        csv_path="example_compounds.csv",
        output_dir="./example_output_basic",
        target_r2=0.7,
        max_iterations=5,
        use_msa_server=False,  # Set to True if MSA server is available
        log_level="INFO"
    )

    # Run the iterator
    results = iterator.run()

    # Print summary
    print("\n" + "="*60)
    print("Results Summary:")
    print("="*60)
    print(f"Converged: {results['converged']}")
    print(f"Final R²: {results['final_r2']:.4f}")
    print(f"Iterations: {results['iterations_run']}/{results['max_iterations']}")
    print(f"Successful predictions: {results['successful_predictions']}/{results['total_compounds']}")


def example_with_msa():
    """
    Example with pre-computed MSA file.
    """
    print("\n" + "="*60)
    print("Example 2: With Pre-computed MSA")
    print("="*60)

    protein_sequence = (
        "STNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIK"
        "KRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEETEIMIVQAKGRGRGRK"
    )

    iterator = Boltz2SARIterator(
        protein_sequence=protein_sequence,
        csv_path="example_compounds.csv",
        output_dir="./example_output_msa",
        target_r2=0.8,
        max_iterations=10,
        use_msa_server=False,
        msa_path="./protein_msa.a3m",  # Provide your MSA file here
        log_level="INFO"
    )

    results = iterator.run()
    return results


def example_custom_parameters():
    """
    Example with custom parameters for optimization.
    """
    print("\n" + "="*60)
    print("Example 3: Custom Parameters")
    print("="*60)

    protein_sequence = (
        "STNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIK"
        "KRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEETEIMIVQAKGRGRGRK"
    )

    iterator = Boltz2SARIterator(
        protein_sequence=protein_sequence,
        csv_path="example_compounds.csv",
        output_dir="./example_output_custom",
        target_r2=0.85,  # Higher R² target
        max_iterations=20,  # More iterations allowed
        use_msa_server=False,
        log_level="DEBUG"  # More detailed logging
    )

    results = iterator.run()
    return results


def example_analyze_results():
    """
    Example of analyzing results after running.
    """
    import pandas as pd
    import json
    from pathlib import Path

    output_dir = Path("./example_output_basic")

    # Load final results
    results_file = output_dir / "final_results.csv"
    if results_file.exists():
        df = pd.read_csv(results_file)
        print("\n" + "="*60)
        print("Compound Results:")
        print("="*60)
        print(df.to_string(index=False))

        # Calculate additional statistics
        df_valid = df.dropna(subset=['Predicted_Affinity'])
        if len(df_valid) > 0:
            print("\nStatistics:")
            print(f"Mean experimental activity: {df_valid['Experimental_Activity'].mean():.2f}")
            print(f"Mean predicted affinity: {df_valid['Predicted_Affinity'].mean():.2f}")
            print(f"Correlation: {df_valid['Experimental_Activity'].corr(df_valid['Predicted_Affinity']):.3f}")

    # Load iteration history
    iterations_file = output_dir / "iteration_history.csv"
    if iterations_file.exists():
        df_iter = pd.read_csv(iterations_file)
        print("\n" + "="*60)
        print("Iteration History:")
        print("="*60)
        print(df_iter.to_string(index=False))

    # Load summary
    summary_file = output_dir / "summary.json"
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            summary = json.load(f)
        print("\n" + "="*60)
        print("Summary:")
        print("="*60)
        for key, value in summary.items():
            print(f"{key}: {value}")


def example_batch_processing():
    """
    Example of processing multiple protein targets.
    """
    print("\n" + "="*60)
    print("Example 4: Batch Processing Multiple Targets")
    print("="*60)

    # Define multiple protein targets
    targets = {
        "target1": "STNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYK",
        # Add more targets as needed
    }

    results_summary = {}

    for target_name, sequence in targets.items():
        print(f"\nProcessing {target_name}...")

        # You would have separate CSV files for each target
        csv_file = f"{target_name}_compounds.csv"

        try:
            iterator = Boltz2SARIterator(
                protein_sequence=sequence,
                csv_path=csv_file,
                output_dir=f"./output_{target_name}",
                target_r2=0.7,
                max_iterations=10,
                use_msa_server=False,
                log_level="INFO"
            )

            results = iterator.run()
            results_summary[target_name] = results

        except Exception as e:
            print(f"Error processing {target_name}: {e}")
            results_summary[target_name] = {"error": str(e)}

    print("\n" + "="*60)
    print("Batch Processing Summary:")
    print("="*60)
    for target, results in results_summary.items():
        if "error" in results:
            print(f"{target}: ERROR - {results['error']}")
        else:
            print(f"{target}: R²={results['final_r2']:.3f}, Converged={results['converged']}")


if __name__ == '__main__':
    # Run basic example
    # NOTE: Make sure Boltz-2 is installed before running:
    # pip install boltz[cuda] -U

    print("Boltz-2 SAR Iterator - Example Usage")
    print("="*60)
    print("This script demonstrates various usage patterns.")
    print("Make sure you have:")
    print("1. Boltz-2 installed (pip install boltz[cuda] -U)")
    print("2. example_compounds.csv file available")
    print("3. Appropriate compute resources (GPU recommended)")
    print("="*60)

    choice = input("\nWhich example would you like to run?\n"
                   "1. Basic usage\n"
                   "2. With MSA (requires MSA file)\n"
                   "3. Custom parameters\n"
                   "4. Analyze existing results\n"
                   "5. Batch processing (requires multiple CSV files)\n"
                   "Enter choice (1-5) or 'q' to quit: ")

    if choice == '1':
        example_basic_usage()
    elif choice == '2':
        example_with_msa()
    elif choice == '3':
        example_custom_parameters()
    elif choice == '4':
        example_analyze_results()
    elif choice == '5':
        example_batch_processing()
    elif choice.lower() == 'q':
        print("Exiting...")
    else:
        print("Invalid choice. Running basic example by default...")
        example_basic_usage()

    print("\nDone!")
