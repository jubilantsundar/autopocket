#!/usr/bin/env python3
"""Quick script to verify correlation calculations for iteration 1"""

import numpy as np
from scipy.stats import spearmanr, pearsonr

# Iteration 1 data (sorted by compound ID)
# IC50 values in nM
experimental = [370, 2670, 7280, 13600, 21200, 25100, 13200]
compound_ids = ['50628047', '50628048', '50628049', '50628050', '50628051', '50628052', '50628053']

# affinity_pred_value from Boltz-2 (log10(IC50) in μM)
predicted = [0.439978688955307, 0.7333598732948303, 0.2827712297439575,
             0.37104663252830505, 0.9437280893325806, 1.120896577835083,
             0.7156575322151184]

print("=" * 60)
print("Iteration 1 Data:")
print("=" * 60)
for i, cid in enumerate(compound_ids):
    print(f"  {cid}: IC50={experimental[i]:>6} nM, Predicted={predicted[i]:.4f}")

print("\n" + "=" * 60)
print("Correlation Metrics:")
print("=" * 60)

# Calculate Pearson R² (as in the code)
exp_array = np.array(experimental)
pred_array = np.array(predicted)

pearson_r, _ = pearsonr(exp_array, pred_array)
r2 = pearson_r ** 2

print(f"  Pearson R: {pearson_r:.4f}")
print(f"  R² (Pearson²): {r2:.4f}")

# Calculate Spearman
spearman_r, _ = spearmanr(exp_array, pred_array)
print(f"  Spearman R: {spearman_r:.4f}")

# Combined metric
combined = max(r2, spearman_r ** 2)
print(f"  Combined metric: {combined:.4f}")

print("\n" + "=" * 60)
print("BEFORE FIX - from log (line 27):")
print("  R² (Pearson²): -1.9691  ❌ WRONG")
print("  Spearman R: 0.5357      ✓ Correct")
print("  Combined metric: 0.2870 ❌ WRONG (using Spearman² since R² was negative)")
print("\n" + "=" * 60)
print("AFTER FIX - Expected values:")
print(f"  R² (Pearson²): {r2:.4f}  ✓ Should be Pearson²")
print(f"  Spearman R: {spearman_r:.4f}      ✓ Already correct")
print(f"  Combined metric: {combined:.4f} ✓ max(R², Spearman²) = max({r2:.4f}, {spearman_r**2:.4f})")
print("=" * 60)
