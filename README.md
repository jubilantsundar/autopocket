# AutoPocket

**Automated binding pocket discovery using AI-powered cofolding**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Powered by Boltz-2](https://img.shields.io/badge/Powered%20by-Boltz--2-green.svg)](https://github.com/jwohlwend/boltz)

> **Transform ligand-based drug design into structure-based design**
> AutoPocket discovers binding pockets that explain your SAR data, enabling structure-based drug design even when pocket information is unknown.

---

## 🎯 Overview

AutoPocket bridges the gap between **ligand-based (LBDD)** and **structure-based drug design (SBDD)** by:

1. **Starting with SAR data** (compounds + activity values)
2. **Discovering binding pockets** using Boltz-2 AI cofolding
3. **Correlating predictions** with experimental activity
4. **Enabling structure-based design** without prior structural knowledge

### Key Innovation

Unlike traditional docking methods that require known binding pockets, **AutoPocket automatically discovers the pocket that best explains your SAR data** through iterative AI-powered cofolding simulations.

---

## ✨ Features

- 🤖 **AI-Powered**: Leverages Boltz-2 generative AI for protein-ligand cofolding
- 🎯 **Blind Pocket Finding**: No prior pocket knowledge required
- 📊 **Multi-Metric Correlation**: Uses R², Spearman, and optional ROC-AUC
- 🔄 **Iterative Refinement**: Converges to optimal pocket predictions
- 📈 **SAR-Driven**: Finds pockets that explain experimental activity
- ⚡ **GPU Accelerated**: Optimized for CUDA-enabled GPUs
- 🌐 **Colab Ready**: Includes Google Colab notebook for cloud execution
- 📁 **Comprehensive Output**: Structures, affinities, correlation metrics, and visualizations

---

## 📦 Installation

### Prerequisites

- Python 3.8 or higher
- CUDA-compatible GPU (recommended, 1000x faster than CPU)
- Git

### Install Dependencies

```bash
# Clone the repository
git clone https://github.com/Jubilantsundar/autopocket.git
cd autopocket

# Install Boltz-2 with CUDA support
pip install boltz[cuda] -U

# Or for CPU-only (slower):
# pip install boltz -U

# Install AutoPocket dependencies
pip install -r requirements.txt
```

---

## 🚀 Quick Start

### 1. Prepare Your Data

Create a CSV file with your SAR data:

```csv
ID,SMILES,Activity
compound_1,Nc1nc(-c2ccc3c(N)n[nH]c3c2)ccn1,370
compound_2,Nc1n[nH]c2cc(-c3ncncc3)ccc12,2670
compound_3,Nc1nccc(-c2ccc3c(N)n[nH]c3c2)c1,7280
```

**Required columns:**
- `SMILES`: Chemical structure in SMILES notation
- `Activity`: Experimental values (e.g., IC50 in nM, pKi, etc.)
- `ID` (optional): Compound identifier

### 2. Run AutoPocket

```bash
python autopocket.py \
  --protein "YOUR_PROTEIN_SEQUENCE" \
  --csv your_compounds.csv \
  --target-r2 0.7 \
  --max-iterations 10
```

### 3. View Results

```
autopocket_output/
├── final_results.csv          # Predictions + experimental values
├── iteration_history.csv      # Convergence tracking
├── summary.json              # Overall statistics
├── predictions/              # Boltz-2 structures and affinities
└── visualizations/           # Correlation plots
```

---

## 📚 Usage Examples

### Basic Usage (Command Line)

```bash
python autopocket.py \
  --protein "MVTPEGNVSLVDESLLVGVTDEDRAV..." \
  --csv compounds.csv \
  --output-dir ./my_results \
  --target-r2 0.7 \
  --max-iterations 10
```

### With Configuration File

```bash
python autopocket.py --config config.json
```

**config.json:**
```json
{
  "protein_sequence": "MVTPEGNVSLVDESLLVGVTDEDRAV...",
  "csv_path": "compounds.csv",
  "output_dir": "./autopocket_output",
  "target_r2": 0.7,
  "max_iterations": 10,
  "protein_chain_id": "A",
  "ligand_chain_id": "L"
}
```

### Python API

```python
from autopocket import AutoPocket

# Initialize
pocket_finder = AutoPocket(
    protein_sequence="MVTPEGNVSLVDESLLVGVTDEDRAV...",
    csv_path="compounds.csv",
    target_r2=0.7,
    max_iterations=10
)

# Run discovery
results = pocket_finder.run()

# Access results
print(f"Final correlation: {results['combined_metric']:.4f}")
print(f"Converged: {results['converged']}")
print(f"Iterations: {results['iterations_run']}")
```

### Advanced: With Constraints

```python
pocket_finder = AutoPocket(
    protein_sequence="...",
    csv_path="compounds.csv",
    pocket_residues=[("A", 107), ("A", 98)],  # Optional hints
    roc_threshold=10000,  # IC50 < 10μM = binder
    log_level="DEBUG"
)
```

---

## 📊 Understanding the Output

### Correlation Metrics

AutoPocket calculates **multiple metrics** and uses the **best one** for convergence:

- **R² (Pearson²)**: Linear correlation
- **Spearman R**: Rank correlation (better for monotonic relationships)
- **Combined**: max(R², Spearman²) - used for convergence
- **ROC-AUC**: Binary classification (optional)

### Affinity Predictions

Boltz-2 outputs:
- `affinity_pred_value`: log₁₀(IC50) in μM units
- `affinity_probability_binary`: Binder probability (0-1)

**Conversion to pIC50:**
```python
pIC50 = 6 - affinity_pred_value
```

---

## 🌐 Google Colab

Run AutoPocket in the cloud with GPU acceleration:

1. Upload `colab_files.zip` to Google Colab
2. Open `AutoPocket_Colab.ipynb`
3. Set runtime to **A100 GPU**
4. Follow the notebook instructions

**Estimated runtime** (A100 GPU):
- 5 compounds, 5 iterations: ~10-15 minutes
- 10 compounds, 10 iterations: ~40-60 minutes

---

## 🔬 How It Works

```
SAR Data → AutoPocket → Binding Pocket → Structure-Based Design
   ↓            ↓              ↓                    ↓
Ligands   AI Cofolding   Correlates with    3D Structures
+ IC50    (Boltz-2)      Experimental       + Affinities
                         Activity
```

**Algorithm:**

1. Load SAR data (SMILES + activity values)
2. For each compound:
   - Generate protein-ligand cofolding structure
   - Extract predicted binding affinity
3. Calculate correlation between predicted and experimental values
4. Iterate until correlation threshold is met
5. Output structures, affinities, and correlation metrics

---

## 📖 Documentation

- [Quick Start Guide](QUICKSTART.md)
- [Configuration Options](docs/configuration.md)
- [API Reference](docs/api.md)
- [Examples](examples/)
- [Troubleshooting](docs/troubleshooting.md)

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

**Note:** Boltz-2 has its own license. Please refer to the [Boltz-2 repository](https://github.com/jwohlwend/boltz) for details.

---

## 📝 Citation

If you use AutoPocket in your research, please cite:

```bibtex
@software{autopocket2024,
  title={AutoPocket: Automated Binding Pocket Discovery for SAR-Driven Drug Design},
  author={Jubilant, Sundar},
  year={2024},
  url={https://github.com/Jubilantsundar/autopocket}
}
```

**Also cite Boltz-2:**
```bibtex
@article{boltz2024,
  title={Accurate Biomolecular Structure Prediction at Scale},
  author={Wohlwend, Jeremy and others},
  journal={bioRxiv},
  year={2024}
}
```

---

## 🙏 Acknowledgments

- **Boltz-2** team for the amazing AI cofolding model
- Inspiration from AlphaFold, RoseTTAFold, and traditional docking methods
- The computational chemistry and drug discovery community

---

## 📧 Support

- **Issues**: [GitHub Issues](https://github.com/Jubilantsundar/autopocket/issues)
- **Discussions**: [GitHub Discussions](https://github.com/Jubilantsundar/autopocket/discussions)
- **Email**: jubilantsundar@gmail.com

---

## ⭐ Star History

If you find AutoPocket useful, please consider giving it a star! ⭐

---

**Made with ❤️ for the drug discovery community**
