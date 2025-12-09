#!/usr/bin/env python3
"""
Visualization script for Boltz-2 SAR Iterator results.

Creates plots showing:
1. Predicted vs Experimental correlation
2. Iteration convergence (R² over iterations)
3. Residual analysis
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
import argparse
from scipy import stats


def plot_correlation(df_results: pd.DataFrame, output_dir: Path):
    """Plot predicted vs experimental activity with correlation metrics."""
    # Filter out None values
    df_plot = df_results.dropna(subset=['Predicted_Affinity'])

    if len(df_plot) < 2:
        print("Not enough data points for correlation plot")
        return

    experimental = df_plot['Experimental_Activity'].values
    predicted = df_plot['Predicted_Affinity'].values

    # Calculate statistics
    r2 = 1 - (np.sum((experimental - predicted) ** 2) /
              np.sum((experimental - np.mean(experimental)) ** 2))
    pearson_r, pearson_p = stats.pearsonr(experimental, predicted)
    spearman_r, spearman_p = stats.spearmanr(experimental, predicted)

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))

    # Scatter plot
    ax.scatter(experimental, predicted, alpha=0.6, s=100, edgecolors='black')

    # Add labels for each point
    for idx, row in df_plot.iterrows():
        ax.annotate(
            row['Compound_ID'],
            (row['Experimental_Activity'], row['Predicted_Affinity']),
            xytext=(5, 5),
            textcoords='offset points',
            fontsize=8,
            alpha=0.7
        )

    # Regression line
    z = np.polyfit(experimental, predicted, 1)
    p = np.poly1d(z)
    x_line = np.linspace(experimental.min(), experimental.max(), 100)
    ax.plot(x_line, p(x_line), "r--", alpha=0.8, linewidth=2, label='Linear fit')

    # Ideal line (y=x)
    min_val = min(experimental.min(), predicted.min())
    max_val = max(experimental.max(), predicted.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--',
            alpha=0.3, linewidth=1, label='Ideal (y=x)')

    # Labels and title
    ax.set_xlabel('Experimental Activity', fontsize=12, fontweight='bold')
    ax.set_ylabel('Predicted Affinity (log10(IC50) μM)', fontsize=12, fontweight='bold')
    ax.set_title('Boltz-2 SAR Correlation Analysis', fontsize=14, fontweight='bold')

    # Add statistics text
    stats_text = f'R² = {r2:.3f}\n'
    stats_text += f'Pearson r = {pearson_r:.3f} (p={pearson_p:.3e})\n'
    stats_text += f'Spearman ρ = {spearman_r:.3f} (p={spearman_p:.3e})\n'
    stats_text += f'N = {len(df_plot)}'

    ax.text(0.05, 0.95, stats_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    output_file = output_dir / 'correlation_plot.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved correlation plot to {output_file}")
    plt.close()


def plot_iteration_history(df_iterations: pd.DataFrame, output_dir: Path, target_r2: float):
    """Plot R² convergence over iterations."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # R² over iterations
    ax1.plot(df_iterations['iteration'], df_iterations['r2'],
             'o-', linewidth=2, markersize=8, label='R²')
    ax1.axhline(y=target_r2, color='r', linestyle='--',
                linewidth=2, alpha=0.7, label=f'Target R² ({target_r2})')
    ax1.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax1.set_ylabel('R²', fontsize=12, fontweight='bold')
    ax1.set_title('R² Convergence', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Success/Failure counts
    width = 0.35
    x = df_iterations['iteration'].values
    ax2.bar(x - width/2, df_iterations['successful'], width,
            label='Successful', alpha=0.8, color='green')
    ax2.bar(x + width/2, df_iterations['failed'], width,
            label='Failed', alpha=0.8, color='red')
    ax2.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax2.set_title('Prediction Success/Failure', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    output_file = output_dir / 'iteration_history.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved iteration history to {output_file}")
    plt.close()


def plot_residuals(df_results: pd.DataFrame, output_dir: Path):
    """Plot residual analysis."""
    df_plot = df_results.dropna(subset=['Predicted_Affinity'])

    if len(df_plot) < 2:
        print("Not enough data points for residual plot")
        return

    experimental = df_plot['Experimental_Activity'].values
    predicted = df_plot['Predicted_Affinity'].values
    residuals = experimental - predicted

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Residual vs Predicted
    ax1.scatter(predicted, residuals, alpha=0.6, s=100, edgecolors='black')
    ax1.axhline(y=0, color='r', linestyle='--', linewidth=2, alpha=0.7)
    ax1.set_xlabel('Predicted Affinity', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Residuals', fontsize=12, fontweight='bold')
    ax1.set_title('Residual Plot', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Residual distribution
    ax2.hist(residuals, bins=min(20, len(residuals)),
             alpha=0.7, edgecolor='black', color='skyblue')
    ax2.axvline(x=0, color='r', linestyle='--', linewidth=2, alpha=0.7)
    ax2.set_xlabel('Residuals', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title('Residual Distribution', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')

    # Add statistics
    stats_text = f'Mean = {np.mean(residuals):.3f}\n'
    stats_text += f'Std Dev = {np.std(residuals):.3f}\n'
    stats_text += f'RMSE = {np.sqrt(np.mean(residuals**2)):.3f}'
    ax2.text(0.95, 0.95, stats_text,
            transform=ax2.transAxes,
            fontsize=10,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    output_file = output_dir / 'residual_analysis.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved residual analysis to {output_file}")
    plt.close()


def create_summary_report(output_dir: Path):
    """Create all visualization plots from Boltz-2 SAR Iterator output."""
    output_dir = Path(output_dir)

    # Load results
    results_file = output_dir / "final_results.csv"
    iterations_file = output_dir / "iteration_history.csv"

    if not results_file.exists():
        print(f"Error: {results_file} not found")
        return

    df_results = pd.read_csv(results_file)
    print(f"Loaded {len(df_results)} compounds from {results_file}")

    # Load summary to get target R²
    summary_file = output_dir / "summary.json"
    target_r2 = 0.7  # default
    if summary_file.exists():
        import json
        with open(summary_file, 'r') as f:
            summary = json.load(f)
            target_r2 = summary.get('target_r2', 0.7)

    # Create plots
    print("\nGenerating visualizations...")
    plot_correlation(df_results, output_dir)
    plot_residuals(df_results, output_dir)

    if iterations_file.exists():
        df_iterations = pd.read_csv(iterations_file)
        plot_iteration_history(df_iterations, output_dir, target_r2)

    print("\nVisualization complete!")


def main():
    parser = argparse.ArgumentParser(
        description='Visualize Boltz-2 SAR Iterator results'
    )
    parser.add_argument(
        'output_dir',
        help='Path to Boltz-2 SAR Iterator output directory'
    )

    args = parser.parse_args()
    create_summary_report(args.output_dir)


if __name__ == '__main__':
    main()
