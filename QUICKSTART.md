# Quick Start Guide - Boltz-2 SAR Iterator

Get started with the Boltz-2 SAR Iterator in 5 minutes!

## Prerequisites

- Python 3.8+
- CUDA-capable GPU (recommended but not required)
- Internet connection (for MSA server, optional)

## Installation

```bash
# 1. Install Boltz-2 with GPU support
pip install boltz[cuda] -U

# Or for CPU-only (slower):
# pip install boltz -U

# 2. Install other dependencies
pip install pandas numpy pyyaml matplotlib seaborn scipy
```

## Quick Test Run

### Step 1: Prepare Your Data

Create a CSV file named `my_compounds.csv`:

```csv
ID,SMILES,Activity
comp1,CC(C)Cc1ccc(cc1)C(C)C(=O)O,6.2
comp2,N[C@@H](Cc1ccc(O)cc1)C(=O)O,7.5
comp3,CC(C)C[C@H](N)C(=O)O,5.8
```

**Required columns:**
- `SMILES`: Chemical structure in SMILES format
- `Activity`: Experimental activity (pIC50, pKi, etc.)
- `ID`: Compound identifier (optional, auto-generated if missing)

### Step 2: Run the Tool

```bash
python boltz2_sar_iterator.py \
  --protein "YOUR_PROTEIN_SEQUENCE_HERE" \
  --csv my_compounds.csv \
  --target-r2 0.7 \
  --max-iterations 5
```

Replace `YOUR_PROTEIN_SEQUENCE_HERE` with your target protein's amino acid sequence.

### Step 3: View Results

Results are saved in `./boltz2_sar_output/`:

```bash
# View final predictions
cat boltz2_sar_output/final_results.csv

# View iteration progress
cat boltz2_sar_output/iteration_history.csv

# View summary
cat boltz2_sar_output/summary.json
```

### Step 4: Visualize Results (Optional)

```bash
python visualize_results.py ./boltz2_sar_output
```

This creates plots:
- `correlation_plot.png`: Predicted vs Experimental
- `iteration_history.png`: R² convergence
- `residual_analysis.png`: Error distribution

## Using the Example Data

We've included example data to test the tool:

```bash
# Using the provided example compounds
python boltz2_sar_iterator.py \
  --protein "STNPPPPETSNPNKPKRQTNQLQYLLRVVLKTLWKHQFAWPFQQPVDAVKLNLPDYYKIIKTPMDMGTIKKRLENNYYWNAQECIQDFNTMFTNCYIYNKPGDDIVLMAEALEKLFLQKINELPTEETEIMIVQAKGRGRGRK" \
  --csv example_compounds.csv \
  --output-dir ./test_run \
  --target-r2 0.7 \
  --max-iterations 3 \
  --log-level INFO
```

## Python API Usage

```python
from boltz2_sar_iterator import Boltz2SARIterator

# Initialize
iterator = Boltz2SARIterator(
    protein_sequence="YOUR_PROTEIN_SEQUENCE",
    csv_path="my_compounds.csv",
    target_r2=0.7,
    max_iterations=5
)

# Run
results = iterator.run()

# Check results
print(f"R²: {results['final_r2']:.3f}")
print(f"Converged: {results['converged']}")
```

## Common Parameters

| Parameter | What it does | Typical values |
|-----------|-------------|----------------|
| `--target-r2` | Stop when R² reaches this value | 0.5 - 0.9 |
| `--max-iterations` | Maximum iterations to run | 5 - 20 |
| `--use-msa-server` | Use online MSA generation | True (default) |
| `--msa-path` | Use pre-computed MSA | Path to .a3m file |
| `--log-level` | Amount of logging detail | INFO, DEBUG |

## Expected Runtime

Typical runtime per compound:
- **With GPU**: 30 seconds - 2 minutes
- **CPU only**: 10 - 30 minutes

For 10 compounds with 5 iterations:
- **With GPU**: ~10-30 minutes
- **CPU only**: ~2-5 hours

## Interpreting Results

### R² Values
- **R² > 0.7**: Good correlation, model is predictive
- **R² 0.5-0.7**: Moderate correlation, some predictive power
- **R² < 0.5**: Poor correlation, check data quality

### Affinity Values
Boltz-2 outputs `log10(IC50)` in μM units:
- Lower values = stronger binding
- Typical range: 0 to 10
- Convert to pIC50: `pIC50 = 6 - affinity_pred_value`

## Troubleshooting

**Problem**: Command not found
```bash
# Make sure boltz is installed
pip show boltz

# Try reinstalling
pip install boltz[cuda] -U --force-reinstall
```

**Problem**: GPU not detected
```bash
# Check CUDA availability
python -c "import torch; print(torch.cuda.is_available())"

# If False, install CPU version or fix CUDA setup
pip install boltz -U  # CPU version
```

**Problem**: CSV parsing error
- Check column names are exactly `SMILES` and `Activity` (case-insensitive OK)
- Ensure no empty cells in required columns
- Verify SMILES strings are valid

**Problem**: All predictions fail
- Check protein sequence is valid (standard amino acids)
- Verify SMILES strings are chemically valid
- Check you have internet (if using MSA server)

## Next Steps

1. **Read the full README.md** for detailed documentation
2. **Run example_usage.py** for more advanced examples
3. **Check the logs** in `boltz2_sar_output/boltz2_sar_iterator.log`
4. **Experiment with parameters** to optimize for your data

## Getting Help

- Check the logs: `boltz2_sar_output/boltz2_sar_iterator.log`
- Refer to Boltz-2 docs: https://github.com/jwohlwend/boltz
- Review example_usage.py for code examples

## Tips for Best Results

1. **Data quality matters**: Clean, well-curated SAR data gives better correlations
2. **Scale matters**: Ensure activity values are on appropriate scale (pIC50, pKi)
3. **Start small**: Test with 5-10 compounds first
4. **Use GPU**: Dramatically faster than CPU
5. **Pre-compute MSA**: Faster for multiple runs with same protein

Happy predicting! 🚀
