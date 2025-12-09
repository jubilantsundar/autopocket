# AutoPocket SAR Workflow - Updated Version

This updated workflow script (`autopocket_sar_workflow_updated.py`) supports all four molecular docking/cofolding tools with proper environment handling.

## Features

- **Boltz-2**: Protein-ligand cofolding with affinity prediction
- **DiffDock**: Molecular docking using separate conda environment
- **Chai-1**: Multi-modal structure prediction
- **ProtoBind-Diff**: Structure-free ligand generation

## Installation Summary

### Tools Installed

1. **Boltz-2** (System Python 3.12)
   ```bash
   cd /content/autopocket/boltz
   pip install -e .[cuda]
   ```

2. **DiffDock** (Conda environment: diffdock, Python 3.9)
   ```bash
   conda env create -f /content/autopocket/DiffDock/environment.yml
   conda activate diffdock
   pip install torch==1.13.1+cu117 --extra-index-url https://download.pytorch.org/whl/cu117
   ```

3. **Chai-1** (System Python 3.12)
   ```bash
   cd /content/autopocket/chai-lab
   pip install -e .
   ```

4. **ProtoBind-Diff** (System Python 3.12)
   ```bash
   cd /content/autopocket/ProtoBind-Diff
   pip install -e .
   ```

## Input Files

### 1. Protein FASTA File
```
>protein_name
FKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVD...
```

### 2. Ligands CSV File
Must contain columns: `compound_id`, `smiles`, `activity`

Example (`example_ligands.csv`):
```csv
compound_id,smiles,activity
erlotinib,Fc1ccc(Nc2ncnc3cc(OCCCN4CCOCC4)c(NC(=O)C=C)cc23)cc1Cl,0.5
gefitinib,COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1,1.2
```

- **activity**: IC50 (nM), Kd (nM), or pIC50 values

## Usage

### Basic Usage (Boltz-2)
```bash
python autopocket_sar_workflow_updated.py \
    protein.fasta \
    ligands.csv \
    --method boltz \
    --correlation 0.7 \
    --output results
```

### Using DiffDock
```bash
python autopocket_sar_workflow_updated.py \
    protein.fasta \
    ligands.csv \
    --method diffdock \
    --correlation 0.7
```

### Using All Methods
```bash
python autopocket_sar_workflow_updated.py \
    protein.fasta \
    ligands.csv \
    --method all \
    --max-iterations 50
```

## Command-Line Options

- `protein`: Path to FASTA file with protein sequence
- `ligands`: Path to CSV file with ligand data
- `--method`: Choose method: `boltz`, `diffdock`, `protobind`, `chai1`, or `all` (default: boltz)
- `--correlation`: Target correlation threshold (default: 0.7)
- `--max-iterations`: Maximum random seed iterations (default: 100)
- `--output`: Output directory (default: sar_workflow_output)
- `--base-dir`: Base installation directory (default: /content/autopocket)

## How It Works

1. **Input Loading**: Loads protein sequence and ligand dataset
2. **Iterative Docking**: Runs docking simulations with different random seeds
3. **Correlation Analysis**: Calculates Pearson/Spearman correlation between:
   - Experimental activity (from CSV)
   - Predicted binding affinity (from simulations)
4. **Convergence**: Stops when target correlation is achieved
5. **Output**: Saves results, plots, and summary

## Output Files

```
sar_workflow_output/
├── final_results.csv           # Best correlation results
├── workflow_summary.json       # Summary statistics
├── correlation_plot.png        # Visualization
├── results_iter*.csv          # Intermediate results
└── boltz_*/                   # Individual simulation outputs
    diffdock_*/
    chai1_*/
    protobind_*/
```

## Key Differences from Original

### 1. Environment Handling
- **DiffDock** runs in separate conda environment (`diffdock`)
- Proper conda activation via shell commands
- All other tools use system Python

### 2. Protein Input
- Uses protein sequence directly (no PDB structure needed)
- Boltz-2: Native FASTA support
- DiffDock: Uses `--protein_sequence` with ESMFold
- Chai-1: FASTA format input
- ProtoBind-Diff: Generates ligands from sequence

### 3. Affinity Prediction
- **Boltz-2**: Uses native affinity prediction module
  - Outputs IC50 predictions and binding probabilities
  - Includes confidence scores (iPTM, ligand_iPTM)
- **DiffDock**: Confidence-based scoring
- **Chai-1**: Aggregate structure quality scores
- **ProtoBind-Diff**: Similarity-based scoring

### 4. Updated Paths
- Configurable base directory
- Automatic tool discovery
- Conda environment activation

## Example: EGFR Kinase SAR Analysis

```bash
# Create protein FASTA
cat > egfr_kinase.fasta << 'EOF'
>EGFR_kinase_domain
FKKIKVLGSGAFGTVYKGLWIPEGEKVKIPVAIKELREATSPKANKEILDEAYVMASVD
NPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYLLNWCVQIAKGMNYLE
DRRLVHRDLAARNVLVKTPQHVKITDFGLAKLLGAEEKEYHAEGGKVPIKWMALESILH
RIYTHQSDVWSYGVTVWELMTFGSKPYDGIPASEISSILEKGERLPQPPICTIDVYMIM
VKCWMIDADSRPKFRELIIEFSKMARDPQRYL
EOF

# Run workflow with Boltz-2
python autopocket_sar_workflow_updated.py \
    egfr_kinase.fasta \
    example_ligands.csv \
    --method boltz \
    --correlation 0.7 \
    --max-iterations 50 \
    --output egfr_sar_results
```

## Troubleshooting

### DiffDock Fails
- Ensure conda environment is activated
- Check PyTorch installation in diffdock env
- Verify CUDA compatibility

### Boltz-2 Affinity Prediction
- Requires `properties` section in YAML
- GPU recommended for speed
- Falls back gracefully on errors

### Memory Issues
- Reduce number of ligands
- Use single method instead of "all"
- Decrease --max-iterations

## Performance Notes

- **Boltz-2**: ~30 seconds per ligand (GPU)
- **DiffDock**: ~2-5 minutes per ligand
- **Chai-1**: ~5-10 minutes per ligand (first run slower due to model downloads)
- **ProtoBind-Diff**: ~1-2 minutes for ligand generation

## Citation

If you use this workflow, please cite:
- Boltz-2: [Paper pending]
- DiffDock: Corso et al., ICLR 2023
- Chai-1: Chai Discovery, bioRxiv 2024
- ProtoBind-Diff: Mistryukova et al., bioRxiv 2025
