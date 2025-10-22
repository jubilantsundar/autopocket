# autopocket
Find a pocket in a protein that explains existing SAR

What This Does
Automatically finds the protein binding pocket that explains your structure-activity relationships (SAR) using cofolding simulations and diffusion-based docking.
No traditional docking required - uses native confidence scores from state-of-the-art methods.

Quick Start
1. Install Dependencies
bash# Install base requirements
pip install numpy scipy pandas matplotlib seaborn rdkit biopython

# Install at least one method (choose based on your needs)
pip install boltz  # Cofolding - good all-around
# OR
pip install chai-lab  # Cofolding - most accurate
# OR
git clone https://github.com/gcorso/DiffDock  # Docking - fastest
cd DiffDock && pip install -r requirements.txt
2. Prepare Input Files
protein.fasta (or just the sequence string):
>MyProtein
MKLLVLLFWLPTLAWQGSEVKPALSELKYPGEYQYRALFQDYQGCYGRGNNSYFLTPEQM...
ligands.csv:
csvsmiles,activity,compound_id
CC(C)Cc1ccc(cc1)C(C)C(=O)O,5.2,ibuprofen
CN1C=NC2=C1C(=O)N(C(=O)N2C)C,230.0,caffeine
COc1ccc2[nH]c3c(c2c1)CCN(C3)C,890.0,melatonin
Required columns:

smiles: Ligand structure in SMILES format
activity: IC50, Kd, or Ki values (lower = better)
compound_id: Unique identifier (optional, will auto-generate)

3. Run Workflow
bash# Basic usage
python sar_workflow.py protein.fasta ligands.csv

# Custom settings
python sar_workflow.py protein.fasta ligands.csv \
    --method boltz \
    --correlation 0.75 \
    --max-iterations 50 \
    --output my_results/
4. Check Results
bash# View summary
cat my_results/workflow_summary.json

# Open correlation plot
open my_results/correlation_plot.png

# Examine detailed results
cat my_results/final_results.csv

Command Line Options
OptionDefaultDescriptionproteinRequiredProtein sequence or FASTA fileligandsRequiredCSV file with SMILES and activity--correlation0.7Target correlation to achieve--max-iterations100Maximum random seed attempts--methodboltzMethod: boltz, alphafold3, chai1, or all--outputsar_workflow_outputOutput directory

Python API Usage
pythonfrom sar_workflow import ProteinLigandWorkflow

# Initialize workflow
workflow = ProteinLigandWorkflow(
    protein_sequence="MKLLVLLF...",  # or "protein.fasta"
    ligand_csv="ligands.csv",
    target_correlation=0.7,
    max_iterations=100,
    output_dir="results/"
)

# Run analysis
results = workflow.run_workflow(method="boltz")

# Access results
print(f"Correlation achieved: {results['correlation']:.3f}")
print(f"Iterations needed: {results['iterations']}")

# Get predictions
df = results['df']
print(df[['compound_id', 'experimental_activity', 'predicted_affinity']])

Understanding the Output
Workflow Summary (JSON)
json{
  "target_correlation": 0.7,
  "achieved_correlation": 0.82,
  "iterations_required": 23,
  "num_ligands": 15,
  "protein_length": 342
}
Results Table (CSV)
compound_idexperimental_activitypredicted_affinityconfidencemethodseedcompound_15.2-85.30.853Boltz-165compound_245.8-68.70.687Boltz-165compound_3230.0-52.10.521Boltz-165
Interpretation:

More negative predicted_affinity = Better predicted binding
Higher confidence = More certain prediction
Same seed = All from same simulation run

Correlation Plot
Shows scatter plot of:

X-axis: Your experimental activity data
Y-axis: Predicted binding scores
Title: Includes Pearson correlation coefficient

Good correlation (r > 0.7) means the model successfully identified the binding pocket.

What Each Method Does
Native Scoring Approach
MethodTypeScore UsedRangeMeaningBoltz-1Cofoldingligand_confidence0-100Confidence in ligand poseDiffDockDockingconfidence0-10+Diffusion model confidenceProtoBind-DiffDockingbinding_energy or confidenceVariableBinding strength estimateChai-1Cofoldingaggregate_score0-1Overall prediction quality
All scores are negated so more negative = better binding (to match IC50/Kd convention).
No traditional docking (Vina, GNINA, etc.) is used - just native scores from the models.
Cofolding vs Docking
Cofolding (Boltz-1, Chai-1):

Input: Protein sequence + Ligand SMILES
Predicts protein structure AND binding simultaneously
Captures protein flexibility
Better for induced-fit binding

Docking (DiffDock, ProtoBind-Diff):

Input: Protein structure + Ligand SMILES
Places ligand into existing structure
Faster (uses pre-generated protein structure)
Better for rigid binding sites


Common Questions
Q: Which method should I use?
For quick screening: Boltz-1 (fastest)
For best accuracy: AlphaFold-3 (requires API key)
For balance: Chai-1 (good speed/accuracy tradeoff)
Q: How many compounds do I need?
Minimum: 8-10 compounds
Recommended: 15-20 compounds
Ideal: 30+ compounds
More compounds = more reliable correlation statistics.
Q: What if correlation is poor?
The workflow automatically tries different random seeds. If after 100 iterations correlation is still < 0.7, consider:

Data quality: Check if IC50 values are from same assay
Binding site: May be allosteric or non-standard
Protein state: May need different conformation
Activity mechanism: May not be driven by binding affinity

Q: Can I use Ki or pIC50 instead of IC50?
Yes! Any activity metric works, as long as:

Consistent units across all compounds
Lower values = better activity (or higher for pIC50/pKi)

For pIC50/pKi, you may want to negate values first since higher pIC50 = better activity.
Q: How long does it take?
Per compound, per iteration:

Boltz-1: ~2-5 minutes (GPU) or ~10-20 minutes (CPU)
AlphaFold-3: ~5-20 minutes (depends on server queue)
Chai-1: ~3-8 minutes (GPU) or ~15-30 minutes (CPU)

Full workflow (15 compounds, ~20 iterations to reach r=0.7):

Boltz-1: ~1-2 hours
AlphaFold-3: ~2-6 hours
Chai-1: ~1.5-3 hours

Q: Do I need a GPU?
Recommended but not required:

GPU: Much faster (10-20x speedup)
CPU: Works but slower

Q: What about AlphaFold-3 API limits?
AlphaFold Server has rate limits. For large-scale screening:

Use Boltz-1 or Chai-1 for initial screening
Validate top hits with AlphaFold-3
Consider local AF3 installation if available


Advanced Features
Run All Methods and Compare
bashpython sar_workflow.py protein.fasta ligands.csv --method all
This runs Boltz-1, AlphaFold-3, AND Chai-1 on all compounds, letting you compare which method performs best for your system.
Custom Correlation Threshold
bash# More stringent
python sar_workflow.py protein.fasta ligands.csv --correlation 0.85

# More lenient
python sar_workflow.py protein.fasta ligands.csv --correlation 0.6
Limit Iterations for Speed
bash# Quick test run
python sar_workflow.py protein.fasta ligands.csv --max-iterations 10

# Exhaustive search
python sar_workflow.py protein.fasta ligands.csv --max-iterations 200

Troubleshooting
Error: "ALPHAFOLD_API_KEY not set"
bash# Set environment variable
export ALPHAFOLD_API_KEY="your_api_key_here"

# Or add to ~/.bashrc or ~/.zshrc
echo 'export ALPHAFOLD_API_KEY="your_key"' >> ~/.bashrc
Error: "CUDA out of memory"
python# In the workflow script, switch to CPU
device="cpu"  # instead of "cuda:0"
Error: "Invalid SMILES"
Check your ligands.csv file:

Remove any invalid characters
Validate SMILES with RDKit
Check for special characters in compound IDs

Poor correlation despite many iterations

Visualize predictions to check if binding site is consistent
Try a different cofolding method
Check experimental data quality
Consider using known binding site constraints


Next Steps
Once you have good correlation (r > 0.7):

Virtual screening: Use the best model to predict binding for new compounds
Lead optimization: Identify which modifications improve binding
Mechanistic insights: Examine predicted structures for key interactions
Experimental validation: Test top predictions in the lab


Support

Documentation: See native_scoring_guide.md for detailed scoring explanation
Integration help: See cofolding_integration.py for API details
Issues: Check GitHub issues or create a new one


Summary
✅ Simple 4-step process: Install → Prepare → Run → Analyze
✅ No docking expertise required
✅ Uses cutting-edge cofolding methods
✅ Automatic optimization via random seed search
✅ Clear visualization and statistics
Happy structure-based drug discovery! 🧬💊
