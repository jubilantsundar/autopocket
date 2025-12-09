#!/bin/bash
# Simple workflow runner to test all tools on random subset of ligands

echo "=========================================="
echo "AutoPocket SAR Workflow"
echo "Comparing protein-ligand binding tools"
echo "=========================================="
echo ""

# Activate conda environment
export PATH="$HOME/miniconda/bin:$PATH"
source $HOME/miniconda/etc/profile.d/conda.sh
conda activate autopocket

# Run workflow with Boltz-2 (using 3 random ligands as a quick test)
echo "Running workflow with Boltz-2..."
python autopocket_sar_workflow_updated.py \
    sequence.txt \
    EGFR-2019-Kd-clean-u-formatted.csv \
    --method boltz \
    --output boltz_output \
    --max-iterations 1

echo ""
echo "=========================================="
echo "Workflow complete! Check boltz_output/ for results"
echo "=========================================="
