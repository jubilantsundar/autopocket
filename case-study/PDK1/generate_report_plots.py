#!/usr/bin/env python3
"""
Generate correlation plots and analysis for PDK1 AutoPocket case study.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr, pearsonr

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

# Load data
output_dir = Path("boltz2_pdk1_sar_output")
results_df = pd.read_csv(output_dir / "final_results.csv")
with open(output_dir / "summary.json") as f:
    summary = json.load(f)

# Calculate metrics
exp_values = results_df['Experimental_Activity'].values
pred_values = results_df['Predicted_Affinity'].values

pearson_r, _ = pearsonr(exp_values, pred_values)
r2 = pearson_r ** 2
spearman_r, _ = spearmanr(exp_values, pred_values)

print(f"Pearson R: {pearson_r:.4f}")
print(f"R²: {r2:.4f}")
print(f"Spearman R: {spearman_r:.4f}")

# Create figure with multiple plots
fig = plt.figure(figsize=(14, 10))

# Plot 1: Scatter plot with regression line
ax1 = plt.subplot(2, 2, 1)
ax1.scatter(exp_values, pred_values, s=100, alpha=0.6, edgecolors='black', linewidth=1.5)

# Add regression line
z = np.polyfit(exp_values, pred_values, 1)
p = np.poly1d(z)
x_line = np.linspace(exp_values.min(), exp_values.max(), 100)
ax1.plot(x_line, p(x_line), "r--", linewidth=2, alpha=0.8, label='Linear fit')

ax1.set_xlabel('Experimental IC50 (nM)', fontsize=12, fontweight='bold')
ax1.set_ylabel('Predicted Affinity (log₁₀ IC50 μM)', fontsize=12, fontweight='bold')
ax1.set_title('PDK1: Predicted vs Experimental Activity', fontsize=14, fontweight='bold')
ax1.legend()

# Add metrics text
textstr = f'R² = {r2:.3f}\nSpearman R = {spearman_r:.3f}\nPearson R = {pearson_r:.3f}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=11,
         verticalalignment='top', bbox=props)

# Plot 2: Rank correlation
ax2 = plt.subplot(2, 2, 2)
ranks_exp = pd.Series(exp_values).rank()
ranks_pred = pd.Series(pred_values).rank()
ax2.scatter(ranks_exp, ranks_pred, s=100, alpha=0.6, edgecolors='black', linewidth=1.5, color='green')

# Add perfect correlation line
ax2.plot([1, 7], [1, 7], 'r--', linewidth=2, alpha=0.8, label='Perfect rank correlation')

ax2.set_xlabel('Experimental Rank', fontsize=12, fontweight='bold')
ax2.set_ylabel('Predicted Rank', fontsize=12, fontweight='bold')
ax2.set_title('Rank Correlation (Spearman)', fontsize=14, fontweight='bold')
ax2.legend()

textstr = f'Spearman R = {spearman_r:.3f}'
ax2.text(0.05, 0.95, textstr, transform=ax2.transAxes, fontsize=11,
         verticalalignment='top', bbox=props)

# Plot 3: Residuals
ax3 = plt.subplot(2, 2, 3)
predicted_from_fit = p(exp_values)
residuals = pred_values - predicted_from_fit
ax3.scatter(exp_values, residuals, s=100, alpha=0.6, edgecolors='black', linewidth=1.5, color='orange')
ax3.axhline(y=0, color='r', linestyle='--', linewidth=2)
ax3.set_xlabel('Experimental IC50 (nM)', fontsize=12, fontweight='bold')
ax3.set_ylabel('Residuals', fontsize=12, fontweight='bold')
ax3.set_title('Residual Plot', fontsize=14, fontweight='bold')

# Plot 4: Bar plot of compounds
ax4 = plt.subplot(2, 2, 4)
compounds = results_df['Compound_ID'].values
x_pos = np.arange(len(compounds))

# Normalize for visualization
exp_norm = (exp_values - exp_values.min()) / (exp_values.max() - exp_values.min())
pred_norm = (pred_values - pred_values.min()) / (pred_values.max() - pred_values.min())

width = 0.35
ax4.bar(x_pos - width/2, exp_norm, width, label='Experimental', alpha=0.8, edgecolor='black')
ax4.bar(x_pos + width/2, pred_norm, width, label='Predicted', alpha=0.8, edgecolor='black')

ax4.set_xlabel('Compound', fontsize=12, fontweight='bold')
ax4.set_ylabel('Normalized Activity', fontsize=12, fontweight='bold')
ax4.set_title('Experimental vs Predicted (Normalized)', fontsize=14, fontweight='bold')
ax4.set_xticks(x_pos)
ax4.set_xticklabels([str(c).replace('50628', '') for c in compounds], rotation=45)
ax4.legend()

plt.tight_layout()
plt.savefig('PDK1_correlation_analysis.png', dpi=300, bbox_inches='tight')
print("\n✓ Saved: PDK1_correlation_analysis.png")

# Create ensemble metrics comparison plot
fig2, axes = plt.subplots(2, 3, figsize=(16, 10))
fig2.suptitle('PDK1: All Metrics Comparison (Iteration 1)', fontsize=16, fontweight='bold')

# Read all metrics from iteration 1
from pathlib import Path
import json

predictions_dir = Path("boltz2_pdk1_sar_output/predictions")
compound_ids = ['50628047', '50628048', '50628049', '50628050', '50628053', '50628051', '50628052']
experimental_ic50 = [370, 2670, 7280, 13600, 13200, 21200, 25100]

metrics_data = {
    'affinity_value': [],
    'probability_binary': [],
    'ensemble_value_mean': [],
    'ensemble_value_std': [],
    'ensemble_prob_mean': [],
    'ensemble_prob_std': []
}

for compound_id in compound_ids:
    pattern = f"boltz_results_iter_1_{compound_id}/predictions/iter_1_{compound_id}/affinity_*.json"
    affinity_files = list(predictions_dir.glob(pattern))

    if affinity_files:
        with open(affinity_files[0]) as f:
            data = json.load(f)

        metrics_data['affinity_value'].append(data.get('affinity_pred_value', np.nan))
        metrics_data['probability_binary'].append(data.get('affinity_probability_binary', np.nan))

        # Calculate ensemble metrics
        ensemble_values = [data.get(f'affinity_pred_value{i}') for i in range(1, 10) if f'affinity_pred_value{i}' in data]
        ensemble_probs = [data.get(f'affinity_probability_binary{i}') for i in range(1, 10) if f'affinity_probability_binary{i}' in data]

        if ensemble_values:
            metrics_data['ensemble_value_mean'].append(np.mean(ensemble_values))
            metrics_data['ensemble_value_std'].append(np.std(ensemble_values))
        else:
            metrics_data['ensemble_value_mean'].append(np.nan)
            metrics_data['ensemble_value_std'].append(np.nan)

        if ensemble_probs:
            metrics_data['ensemble_prob_mean'].append(np.mean(ensemble_probs))
            metrics_data['ensemble_prob_std'].append(np.std(ensemble_probs))
        else:
            metrics_data['ensemble_prob_mean'].append(np.nan)
            metrics_data['ensemble_prob_std'].append(np.nan)

# Plot each metric
metric_names = ['affinity_value', 'probability_binary', 'ensemble_value_mean',
                'ensemble_value_std', 'ensemble_prob_mean', 'ensemble_prob_std']
metric_titles = ['Affinity Value', 'Probability Binary', 'Ensemble Value Mean',
                 'Ensemble Value Std', 'Ensemble Prob Mean', 'Ensemble Prob Std ⭐']

for idx, (metric_name, metric_title) in enumerate(zip(metric_names, metric_titles)):
    ax = axes[idx // 3, idx % 3]

    metric_vals = np.array(metrics_data[metric_name])
    valid_mask = ~np.isnan(metric_vals)

    if valid_mask.sum() > 1:
        exp_subset = np.array(experimental_ic50)[valid_mask]
        metric_subset = metric_vals[valid_mask]

        # Calculate correlations
        try:
            r2_metric = pearsonr(exp_subset, metric_subset)[0] ** 2
            spearman_metric = spearmanr(exp_subset, metric_subset)[0]
        except:
            r2_metric = 0
            spearman_metric = 0

        ax.scatter(exp_subset, metric_subset, s=100, alpha=0.6, edgecolors='black', linewidth=1.5)

        # Add regression line
        if len(exp_subset) > 1:
            z = np.polyfit(exp_subset, metric_subset, 1)
            p = np.poly1d(z)
            x_line = np.linspace(exp_subset.min(), exp_subset.max(), 100)
            ax.plot(x_line, p(x_line), "r--", linewidth=2, alpha=0.8)

        ax.set_xlabel('Experimental IC50 (nM)', fontsize=10)
        ax.set_ylabel(metric_title, fontsize=10)
        ax.set_title(metric_title, fontsize=12, fontweight='bold')

        # Add metrics
        textstr = f'R² = {r2_metric:.3f}\nSpearman = {spearman_metric:.3f}'
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', bbox=props)

plt.tight_layout()
plt.savefig('PDK1_all_metrics_comparison.png', dpi=300, bbox_inches='tight')
print("✓ Saved: PDK1_all_metrics_comparison.png")

print("\n" + "="*70)
print("Report generation complete!")
print("="*70)
