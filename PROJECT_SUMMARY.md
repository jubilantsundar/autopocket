# Boltz-2 SAR Iterator - Project Summary

## Overview

This project provides a comprehensive iterative tool for running Boltz-2 protein-ligand cofolding simulations with Structure-Activity Relationship (SAR) data. The tool automatically correlates predicted binding affinities with experimental activity values until a target R² threshold is achieved.

## Files Created

### Main Tool
- **`boltz2_sar_iterator.py`** - Main iterative tool with full functionality
  - Supports command-line and Python API usage
  - Configurable via command-line arguments or config files
  - Includes chain IDs, pocket constraints, contact constraints, and templates
  - Comprehensive logging and error handling
  - Generates detailed output reports

### Configuration Files
- **`config_example.json`** - Example JSON configuration file
- **`config_example.yaml`** - Example YAML configuration file (with comments)

### Example Data
- **`example_compounds.csv`** - Sample SAR data with 10 compounds
  - Contains SMILES structures and activity values
  - Ready to use for testing

### Utilities
- **`visualize_results.py`** - Results visualization script
  - Creates correlation plots
  - Plots iteration convergence
  - Generates residual analysis
  - Publication-quality figures

- **`example_usage.py`** - Interactive examples demonstrating:
  - Basic usage
  - Usage with pre-computed MSA
  - Custom parameters
  - Result analysis
  - Batch processing multiple targets

### Documentation
- **`README.md`** - Comprehensive documentation
  - Installation instructions
  - Detailed usage examples
  - Input/output file formats
  - Troubleshooting guide
  - Advanced features

- **`QUICKSTART.md`** - Quick start guide
  - 5-minute setup
  - Simple examples
  - Common parameters
  - Tips for best results

- **`requirements.txt`** - Python dependencies
  - Boltz-2 and all required packages
  - Optional visualization packages

- **`PROJECT_SUMMARY.md`** - This file

## Key Features

### 1. Flexible Input Methods
- Command-line arguments for all parameters
- JSON/YAML configuration files
- Python API for programmatic use

### 2. Advanced Boltz-2 Integration
- Custom chain IDs for protein and ligand
- Pocket residue constraints
- Contact residue constraints with distance thresholds
- Template structure support
- MSA server or pre-computed MSA files

### 3. Iterative Correlation
- Runs until target R² is achieved or max iterations reached
- Tracks progress across iterations
- Comprehensive logging of all operations

### 4. Rich Output
- Final predictions CSV with all compounds
- Iteration history tracking R² convergence
- Summary JSON with overall statistics
- Detailed logs for debugging
- Optional visualization plots

## Usage Examples

### Command-Line (Basic)
```bash
python boltz2_sar_iterator.py \
  --protein "MVTPEGNVSL..." \
  --csv compounds.csv \
  --target-r2 0.7 \
  --max-iterations 10
```

### Command-Line (Advanced)
```bash
python boltz2_sar_iterator.py \
  --protein "MVTPEGNVSL..." \
  --csv compounds.csv \
  --protein-chain A \
  --ligand-chain L \
  --pocket-residues "A:107,A:98" \
  --contact-residues "A:20-B:27:5.0,A:21-B:31:5.0" \
  --template template.cif \
  --msa-path protein_msa.a3m \
  --no-msa-server \
  --target-r2 0.8 \
  --max-iterations 15
```

### Using Config File
```bash
python boltz2_sar_iterator.py --config config_example.yaml
```

### Python API
```python
from boltz2_sar_iterator import Boltz2SARIterator

iterator = Boltz2SARIterator(
    protein_sequence="MVTPEGNVSL...",
    csv_path="compounds.csv",
    protein_chain_id="A",
    ligand_chain_id="L",
    pocket_residues=[("A", 107), ("A", 98)],
    contact_residues=[(("A", 20), ("B", 27), 5.0)],
    template_file="template.cif",
    target_r2=0.7,
    max_iterations=10
)

results = iterator.run()
print(f"Final R²: {results['final_r2']:.3f}")
```

### Visualize Results
```bash
python visualize_results.py ./boltz2_sar_output
```

## Configuration Options

### Required Parameters
- `protein_sequence`: Protein amino acid sequence
- `csv_path`: Path to SAR data CSV file

### CSV Format
Required columns:
- `SMILES`: Chemical structure in SMILES notation
- `Activity`: Experimental activity values
- `ID` (optional): Compound identifier

### Optional Parameters
- `output_dir`: Output directory (default: ./boltz2_sar_output)
- `target_r2`: Target R² for convergence (default: 0.7)
- `max_iterations`: Maximum iterations (default: 10)
- `use_msa_server`: Use MSA server (default: true)
- `msa_path`: Pre-computed MSA file path
- `protein_chain_id`: Protein chain ID (default: A)
- `ligand_chain_id`: Ligand chain ID (default: L)
- `pocket_residues`: Binding pocket residues
- `contact_residues`: Contact constraints
- `template_file`: Template CIF file
- `template_force`: Force template usage (default: true)
- `template_threshold`: Template threshold (default: 2)
- `log_level`: Logging detail (default: INFO)

## Output Structure

```
boltz2_sar_output/
├── final_results.csv          # Predictions for all compounds
├── iteration_history.csv      # R² and stats per iteration
├── summary.json              # Overall summary
├── boltz2_sar_iterator.log   # Detailed execution log
├── yaml_inputs/              # Generated Boltz-2 input files
│   └── iter_1_compound_1.yaml
├── predictions/              # Boltz-2 output directories
│   └── iter_1_compound_1/
│       ├── affinity_*.json
│       ├── confidence_*.json
│       └── *.cif
└── [optional plots from visualize_results.py]
    ├── correlation_plot.png
    ├── iteration_history.png
    └── residual_analysis.png
```

## Quick Start

1. **Install dependencies:**
   ```bash
   pip install boltz[cuda] -U
   pip install pandas numpy pyyaml matplotlib seaborn scipy
   ```

2. **Prepare your CSV file with SMILES and Activity columns**

3. **Run the tool:**
   ```bash
   python boltz2_sar_iterator.py \
     --protein "YOUR_SEQUENCE" \
     --csv your_compounds.csv
   ```

4. **View results:**
   ```bash
   cat boltz2_sar_output/summary.json
   python visualize_results.py ./boltz2_sar_output
   ```

## Performance Notes

- **GPU recommended**: ~1000x faster than CPU
- **Typical runtime per compound**: 30s-2min (GPU), 10-30min (CPU)
- **For 10 compounds, 5 iterations**: ~10-30min (GPU), 2-5hrs (CPU)

## Advanced Features

- **Multiple protein targets**: Process different proteins with same compounds
- **Batch processing**: Programmatically run multiple experiments
- **Custom correlation metrics**: Modify code to use Spearman, Kendall, etc.
- **Constraint customization**: Fine-tune binding pocket and contact definitions
- **Template-guided predictions**: Use known structures to guide predictions

## Troubleshooting

See `README.md` for detailed troubleshooting guide covering:
- Installation issues
- CSV parsing errors
- GPU/CUDA problems
- MSA server connection issues
- Low correlation results

## Citation

If using this tool, please cite Boltz-2:
- GitHub: https://github.com/jwohlwend/boltz

## License

This tool is provided for research purposes. Refer to Boltz-2's license for details.

## Development

All scripts are modular and well-documented:
- Clear class structure in `Boltz2SARIterator`
- Comprehensive docstrings
- Type hints throughout
- Logging at multiple levels
- Error handling and validation

## Getting Help

1. Check `QUICKSTART.md` for common tasks
2. Read `README.md` for detailed documentation
3. Review `example_usage.py` for code examples
4. Check logs in `boltz2_sar_output/boltz2_sar_iterator.log`
5. Refer to Boltz-2 documentation

## Future Enhancements

Potential additions:
- Parallel processing of compounds
- Support for multiple protein chains
- Integration with molecular property calculators
- Machine learning-based result refinement
- Web interface
- Docker containerization

---

**Created:** 2025-12-08
**Status:** Production-ready
**Python:** 3.8+
**Dependencies:** See requirements.txt
