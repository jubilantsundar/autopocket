# Google Colab Setup Guide - Boltz-2 SAR Iterator

## Quick Start (5 steps)

### 1. Upload to Google Colab

**Option A: Direct Upload**
1. Go to https://colab.research.google.com
2. File → Upload notebook
3. Upload `Boltz2_SAR_Iterator_Colab.ipynb`

**Option B: From Google Drive**
1. Upload `Boltz2_SAR_Iterator_Colab.ipynb` to your Google Drive
2. Double-click to open in Colab

### 2. Set A100 GPU Runtime

1. Runtime → Change runtime type
2. Hardware accelerator: **GPU**
3. GPU type: **A100** (recommended) or T4
4. Save

### 3. Prepare Files to Upload

You'll need these files from your local machine:
- ✓ `boltz2_sar_iterator.py` (required)
- ✓ Your SAR data CSV file (required)
- ✓ `visualize_results.py` (optional, for plots)
- ✓ MSA file .a3m (optional)
- ✓ Template .pdb or .cif (optional)

### 4. Run the Notebook

Execute cells in order:
1. **Cell 1**: Check GPU (verify A100 is active)
2. **Cell 2-3**: Install dependencies (~5-10 min first time)
3. **Cell 4-5**: Upload files (boltz2_sar_iterator.py and your CSV)
4. **Cell 6**: Configure parameters (edit your protein sequence)
5. **Cell 7**: Run iterator (main execution)
6. **Cell 8-9**: View results
7. **Cell 10**: Download results

### 5. Download Results

The notebook will create a ZIP file with all results. Download it before the session ends.

---

## Detailed Setup Instructions

### Files to Upload (What & When)

#### Required Files (Upload in Cell 3)
```
boltz2_sar_iterator.py  ← The main tool
your_data.csv          ← Your SAR data (SMILES + Activity columns)
```

#### Optional Files (Upload in Cell 4)
```
protein.a3m            ← Pre-computed MSA (if not using MSA server)
template.pdb           ← Template structure (automatically converts to CIF)
visualize_results.py   ← For creating plots (upload before Cell 9)
```

### CSV File Format

Your CSV must have these columns:
```csv
ID,SMILES,Activity
compound_1,CC(C)Cc1ccc(cc1)C(C)C(=O)O,6.2
compound_2,N[C@@H](Cc1ccc(O)cc1)C(=O)O,7.5
compound_3,CC(C)C[C@H](N)C(=O)O,5.8
```

**Required columns:**
- `SMILES`: Chemical structure in SMILES format
- `Activity`: Experimental activity values (pIC50, pKi, etc.)

**Optional columns:**
- `ID`: Compound identifier (auto-generated if missing)

### Configuration (Cell 6)

Edit these variables in Cell 6:

```python
# Required
PROTEIN_SEQUENCE = "YOUR_PROTEIN_SEQUENCE_HERE"

# Adjust as needed
TARGET_R2 = 0.7           # Target correlation
MAX_ITERATIONS = 10       # Maximum iterations
PROTEIN_CHAIN = "A"       # Protein chain ID
LIGAND_CHAIN = "L"        # Ligand chain ID

# Optional constraints
POCKET_RESIDUES = [("A", 107), ("A", 98)]  # Binding pocket
CONTACT_RESIDUES = [                        # Contact constraints
    (("A", 20), ("B", 27), 5.0)
]
```

---

## Runtime Estimates (A100 GPU)

| Compounds | Iterations | Estimated Time |
|-----------|-----------|----------------|
| 5         | 5         | 10-15 min      |
| 10        | 5         | 20-30 min      |
| 10        | 10        | 40-60 min      |
| 20        | 10        | 1.5-2 hours    |
| 50        | 10        | 4-6 hours      |

**Note:** T4 GPU is ~3x slower than A100

---

## Colab Runtime Limits

### Free Tier
- **Runtime**: ~12 hours maximum
- **GPU Access**: Limited hours per week
- **Memory**: 12-16 GB RAM
- **Recommendation**: Process ≤20 compounds at a time

### Colab Pro ($10/month)
- **Runtime**: Up to 24 hours
- **Priority GPU**: Better A100 availability
- **Memory**: Up to 52 GB RAM
- **Recommendation**: Can handle 50+ compounds

### Colab Pro+ ($50/month)
- **Runtime**: Longer sessions
- **Background execution**: Runs even when browser closed
- **Priority hardware**: Best GPU access

---

## Important Colab Tips

### 1. Prevent Disconnection
- Keep the Colab tab active
- Don't close the browser
- Use Colab Pro for background execution
- Consider running overnight

### 2. Save Results Frequently
```python
# Download intermediate results (run this cell periodically)
from google.colab import files
files.download("/content/boltz2_output/summary.json")
files.download("/content/boltz2_output/final_results.csv")
```

### 3. Monitor Progress
```python
# Check GPU usage
!nvidia-smi

# View logs in real-time
!tail -f /content/boltz2_output/boltz2_sar_iterator.log

# Check disk space
!df -h
```

### 4. Resume After Disconnect
If your session disconnects:
1. Reconnect to runtime
2. Upload `boltz2_sar_iterator.py` again
3. Run the download cells to retrieve any saved results
4. Check `/content/boltz2_output/` for partial results

### 5. Memory Management
```python
# Clear output to free memory
from IPython.display import clear_output
clear_output()

# Restart runtime if needed
Runtime → Restart runtime
```

---

## Output Files

After completion, you'll have:

```
boltz2_output/
├── final_results.csv          # Predictions for all compounds
├── iteration_history.csv      # R² convergence
├── summary.json              # Overall statistics
├── boltz2_sar_iterator.log   # Detailed logs
├── correlation_plot.png      # Predicted vs Experimental
├── iteration_history.png     # Convergence plot
├── residual_analysis.png     # Error distribution
└── predictions/              # Detailed Boltz-2 outputs
    └── iter_1_compound_1/
        ├── affinity_*.json
        ├── confidence_*.json
        └── *.cif
```

**Download the ZIP** (Cell 10) to get all files at once.

---

## Troubleshooting

### "No GPU available"
- Runtime → Change runtime type → Select GPU
- Try T4 if A100 unavailable
- Free tier may have limited GPU quota

### "Out of Memory"
```python
# Reduce number of compounds
# Process in smaller batches
# Use simpler protein (<500 aa)
# Clear outputs: Edit → Clear all outputs
```

### "Session disconnected"
```python
# Results are auto-saved to /content/boltz2_output/
# Reconnect and run download cells to retrieve data
# Consider Colab Pro for longer sessions
```

### "MSA server connection failed"
```python
# Option 1: Upload pre-computed MSA file
# Set USE_MSA_SERVER = False
# Upload your .a3m file

# Option 2: Retry with MSA server
# Network issues usually resolve quickly
```

### "Installation failed"
```python
# Restart runtime: Runtime → Restart runtime
# Run installation cells again
# Check internet connection
```

---

## Performance Optimization

### For Faster Execution
1. **Use A100 GPU** (not T4)
2. **Pre-compute MSA** offline, upload .a3m file
3. **Reduce iterations** for initial testing (e.g., 3-5)
4. **Process fewer compounds** per run
5. **Use template structures** when available

### For Large Datasets (50+ compounds)
1. **Split into batches** of 20-30 compounds
2. **Use Colab Pro** for longer runtime
3. **Download results** after each batch
4. **Combine results** offline later

---

## Example: Quick Test Run

Try this for a fast test (should complete in ~5-10 minutes):

```python
# Cell 6 Configuration
PROTEIN_SEQUENCE = "MVTPEGNVSLVDESLLVGVTDEDRAV"  # Short test sequence
TARGET_R2 = 0.5                                   # Lower threshold
MAX_ITERATIONS = 3                                # Fewer iterations

# Test with 3-5 compounds only
```

This validates your setup without using too much GPU time.

---

## Cost Considerations

### Free Tier
- **Cost**: $0
- **Best for**: Testing, small datasets (<10 compounds)
- **Limitations**: GPU quota, session timeouts

### Colab Pro ($10/month)
- **Cost**: $10/month
- **Best for**: Regular use, 10-50 compounds per run
- **Benefits**: Longer runtime, better GPU access

### Colab Pro+ ($50/month)
- **Cost**: $50/month
- **Best for**: Heavy users, large datasets (50+ compounds)
- **Benefits**: Background execution, priority hardware

### Alternative: Local GPU
If you have local GPU (e.g., RTX 4090, A100):
- Use the standalone Python scripts
- No runtime limits
- No upload/download needed
- One-time hardware cost

---

## Support & Resources

### Check Status
```python
# View current configuration
!python boltz2_sar_iterator.py --help

# Check Boltz installation
!boltz --version

# View full logs
!cat /content/boltz2_output/boltz2_sar_iterator.log
```

### Get Help
- **Logs**: Check `boltz2_sar_iterator.log` for errors
- **Boltz-2 Docs**: https://github.com/jwohlwend/boltz
- **Colab FAQ**: https://research.google.com/colaboratory/faq.html

---

## Checklist Before Running

- [ ] A100 GPU runtime selected
- [ ] `boltz2_sar_iterator.py` uploaded
- [ ] CSV file uploaded (with SMILES and Activity columns)
- [ ] Protein sequence configured in Cell 6
- [ ] Parameters adjusted (TARGET_R2, MAX_ITERATIONS)
- [ ] Optional files uploaded (MSA, template, visualize_results.py)
- [ ] Browser tab will stay open during run
- [ ] Aware of estimated runtime for your dataset

---

**Ready to start?** Open the notebook and begin with Cell 1! 🚀
