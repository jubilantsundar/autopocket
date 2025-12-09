#!/usr/bin/env python3
"""Test if all required packages are available"""

print("Testing package imports...\n")

# Core packages
try:
    import numpy as np
    print(f"✓ NumPy: {np.__version__}")
except Exception as e:
    print(f"✗ NumPy: {e}")

try:
    import pandas as pd
    print(f"✓ Pandas: {pd.__version__}")
except Exception as e:
    print(f"✗ Pandas: {e}")

try:
    import scipy
    print(f"✓ SciPy: {scipy.__version__}")
except Exception as e:
    print(f"✗ SciPy: {e}")

# Visualization
try:
    import matplotlib
    print(f"✓ Matplotlib: {matplotlib.__version__}")
except Exception as e:
    print(f"✗ Matplotlib: {e}")

try:
    import seaborn
    print(f"✓ Seaborn: {seaborn.__version__}")
except Exception as e:
    print(f"✗ Seaborn: {e}")

# ML frameworks
try:
    import torch
    print(f"✓ PyTorch: {torch.__version__}")
    print(f"  CUDA available: {torch.cuda.is_available()}")
except Exception as e:
    print(f"✗ PyTorch: {e}")

# Chemistry
try:
    from rdkit import Chem
    print(f"✓ RDKit: Available")
except Exception as e:
    print(f"✗ RDKit: {e}")

# Biology
try:
    import Bio
    print(f"✓ BioPython: {Bio.__version__}")
except Exception as e:
    print(f"✗ BioPython: {e}")

# Docking tools
print("\nDocking tools:")
try:
    import boltz
    print(f"✓ Boltz: Available")
except Exception as e:
    print(f"✗ Boltz: Not installed ({e})")

try:
    from chai_lab import chai1
    print(f"✓ Chai-1: Available")
except Exception as e:
    print(f"✗ Chai-1: Not installed ({e})")

try:
    import transformers
    print(f"✓ Transformers: {transformers.__version__}")
except Exception as e:
    print(f"✗ Transformers: {e}")

try:
    from fair_esm import pretrained
    print(f"✓ ESM: Available")
except Exception as e:
    print(f"✗ ESM: Not installed ({e})")

print("\n" + "="*50)
print("Test complete!")
