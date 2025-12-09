================================================================================
Boltz-2 SAR Iterator - Colab Upload Package
================================================================================

This folder contains everything you need to run Boltz-2 SAR Iterator on
Google Colab with A100 GPU.

FILES IN THIS PACKAGE:
======================

1. Boltz2_SAR_Iterator_Colab.ipynb  ← Upload this to Google Colab
2. boltz2_sar_iterator.py           ← Upload to Colab (Cell 3)
3. visualize_results.py             ← Upload to Colab (Cell 9, optional)
4. example_compounds.csv            ← Sample data (replace with yours)
5. COLAB_SETUP.md                   ← Read this for detailed instructions
6. README_FIRST.txt                 ← This file

QUICK START:
============

Step 1: Upload Notebook to Colab
   - Go to https://colab.research.google.com
   - File → Upload notebook
   - Select: Boltz2_SAR_Iterator_Colab.ipynb

Step 2: Set A100 GPU Runtime
   - Runtime → Change runtime type
   - Hardware accelerator: GPU
   - GPU type: A100
   - Save

Step 3: Prepare Your Data
   - Replace 'example_compounds.csv' with your data
   - Your CSV must have: SMILES and Activity columns
   - Optional: Add MSA file (.a3m) and template (.pdb or .cif)

Step 4: Run the Notebook
   - Execute cells 1-3: Install dependencies
   - Cell 4: Upload boltz2_sar_iterator.py
   - Cell 5: Upload your CSV file
   - Cell 6: Configure your protein sequence
   - Cell 7: Run the iterator
   - Cell 8-10: View and download results

YOUR DATA:
==========

Before uploading to Colab, prepare:
   [ ] Your SAR data CSV (SMILES + Activity columns)
   [ ] Your protein sequence (amino acid sequence)
   [ ] Optional: MSA file (.a3m)
   [ ] Optional: Template structure (.pdb or .cif)

ESTIMATED TIME (A100 GPU):
==========================

5 compounds, 5 iterations:   ~10-15 minutes
10 compounds, 5 iterations:  ~20-30 minutes
10 compounds, 10 iterations: ~40-60 minutes
20 compounds, 10 iterations: ~1.5-2 hours

NEED HELP?
==========

1. Read COLAB_SETUP.md for detailed instructions
2. Check the notebook - it has inline documentation
3. Monitor progress in the notebook output
4. Download results before session ends

TIPS:
=====

✓ Keep browser tab open during execution
✓ Use Colab Pro for longer runtime (>12 hours)
✓ Download results frequently
✓ Start with small test (3-5 compounds, 3 iterations)
✓ A100 is 3x faster than T4

================================================================================
Ready? Upload Boltz2_SAR_Iterator_Colab.ipynb to Colab and follow the steps!
================================================================================
