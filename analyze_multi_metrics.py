#!/usr/bin/env python3
"""
Analyze existing Boltz-2 predictions using multiple metrics.
This script reads affinity JSON files and calculates correlations for all available metrics.
"""

import json
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr, pearsonr
import pandas as pd

# Experimental IC50 values (in nM)
experimental_data = {
    '50628047': 370,
    '50628048': 2670,
    '50628049': 7280,
    '50628050': 13600,
    '50628053': 13200,
    '50628051': 21200,
    '50628052': 25100
}

def extract_all_metrics(affinity_file):
    """Extract all available metrics from an affinity JSON file."""
    with open(affinity_file, 'r') as f:
        data = json.load(f)

    metrics = {}

    # Main predictions
    if 'affinity_pred_value' in data:
        metrics['affinity_value'] = float(data['affinity_pred_value'])
    if 'affinity_probability_binary' in data:
        metrics['probability_binary'] = float(data['affinity_probability_binary'])

    # Ensemble predictions
    ensemble_values = []
    ensemble_probs = []

    for i in range(1, 10):
        value_key = f'affinity_pred_value{i}'
        prob_key = f'affinity_probability_binary{i}'

        if value_key in data:
            ensemble_values.append(float(data[value_key]))
        if prob_key in data:
            ensemble_probs.append(float(data[prob_key]))

    if ensemble_values:
        metrics['ensemble_value_mean'] = np.mean(ensemble_values)
        metrics['ensemble_value_std'] = np.std(ensemble_values)
    if ensemble_probs:
        metrics['ensemble_prob_mean'] = np.mean(ensemble_probs)
        metrics['ensemble_prob_std'] = np.std(ensemble_probs)

    return metrics

def calculate_correlation(exp_values, pred_values):
    """Calculate R² and Spearman for a set of predictions."""
    exp_array = np.array(exp_values)
    pred_array = np.array(pred_values)

    # Pearson R²
    pearson_r, _ = pearsonr(exp_array, pred_array)
    r2 = pearson_r ** 2

    # Spearman
    spearman_r, _ = spearmanr(exp_array, pred_array)

    # Combined
    combined = max(r2, spearman_r ** 2)

    return {
        'r2': r2,
        'pearson_r': pearson_r,
        'spearman_r': spearman_r,
        'combined': combined
    }

def analyze_iteration(predictions_dir, iteration):
    """Analyze all metrics for a specific iteration."""
    print(f"\n{'='*70}")
    print(f"Iteration {iteration} - Multi-Metric Analysis")
    print(f"{'='*70}")

    # Collect all metrics for all compounds
    compound_metrics = {}

    for compound_id in experimental_data.keys():
        # Find affinity file
        pattern = f"boltz_results_iter_{iteration}_{compound_id}/predictions/iter_{iteration}_{compound_id}/affinity_*.json"
        affinity_files = list(predictions_dir.glob(pattern))

        if affinity_files:
            metrics = extract_all_metrics(affinity_files[0])
            compound_metrics[compound_id] = metrics

    if not compound_metrics:
        print("No prediction files found for this iteration")
        return None

    print(f"Found predictions for {len(compound_metrics)}/{len(experimental_data)} compounds\n")

    # Get all unique metric names
    all_metric_names = set()
    for metrics in compound_metrics.values():
        all_metric_names.update(metrics.keys())

    # Calculate correlations for each metric
    results = {}

    for metric_name in sorted(all_metric_names):
        exp_values = []
        pred_values = []

        for compound_id, ic50 in experimental_data.items():
            if compound_id in compound_metrics and metric_name in compound_metrics[compound_id]:
                exp_values.append(ic50)
                pred_values.append(compound_metrics[compound_id][metric_name])

        if len(pred_values) >= 2:
            corr = calculate_correlation(exp_values, pred_values)
            results[metric_name] = corr

            print(f"{metric_name}:")
            print(f"  R² = {corr['r2']:.4f}")
            print(f"  Pearson R = {corr['pearson_r']:.4f}")
            print(f"  Spearman R = {corr['spearman_r']:.4f}")
            print(f"  Combined = {corr['combined']:.4f}")
            print()

    # Find best metric
    if results:
        best_metric = max(results.items(), key=lambda x: x[1]['combined'])
        print(f"✓ BEST METRIC: {best_metric[0]} (Combined = {best_metric[1]['combined']:.4f})")
        print(f"  R² = {best_metric[1]['r2']:.4f}, Spearman = {best_metric[1]['spearman_r']:.4f}")

    return results

if __name__ == "__main__":
    predictions_dir = Path("/Users/sundar/work/my_tools/autopocket/case-study/PDK1/boltz2_pdk1_sar_output/predictions")

    # Analyze available iterations
    available_iterations = set()
    for path in predictions_dir.glob("boltz_results_iter_*"):
        iter_num = int(path.name.split('_')[3])
        available_iterations.add(iter_num)

    print(f"\nFound data for iterations: {sorted(available_iterations)}")

    # Analyze each iteration
    all_results = {}
    for iteration in sorted(available_iterations):
        results = analyze_iteration(predictions_dir, iteration)
        if results:
            all_results[iteration] = results

    # Summary comparison
    print(f"\n{'='*70}")
    print("SUMMARY: Best Metric Per Iteration")
    print(f"{'='*70}")

    for iteration, results in sorted(all_results.items()):
        if results:
            best_metric = max(results.items(), key=lambda x: x[1]['combined'])
            print(f"Iteration {iteration}: {best_metric[0]} (Combined = {best_metric[1]['combined']:.4f})")
