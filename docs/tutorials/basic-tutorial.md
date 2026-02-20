# Basic PRISM Tutorial

!!! example "Quick Start"
    ```bash
    prism protein.pdb ligand.mol2 -o t4l_benzene
    cd t4l_benzene/GMX_PROLIG_MD && bash localrun.sh
    ```

## Overview

This tutorial will guide you through the complete workflow of building a protein-ligand MD system using PRISM, from input preparation to analyzing simulation results.

**Estimated Time:** 30-60 minutes (+ simulation time)

**What You'll Learn:**

- Building an MD-ready system with PRISM
- Running equilibration and production simulations
- Analyzing protein-ligand interactions
- Visualizing contacts and binding modes

## Prerequisites

- PRISM installed ([Installation Guide](../getting-started/installation.md))
- GROMACS 2020+ installed
- Basic command-line knowledge
- ~10 GB free disk space

## Tutorial System

We'll study **T4 Lysozyme L99A/M102Q** with a **benzene** ligand (PDB: [5JWT](https://www.rcsb.org/structure/5JWT)) - a well-characterized model system for protein-ligand binding.

**Why this system?**

- Small protein (164 residues, fast simulations)
- Well-studied hydrophobic binding pocket
- Experimental binding data available for validation
- Benzene is a minimal ligand (12 atoms), ideal for learning

## Step 1: Download Input Files

Pre-processed input files are provided. The protein has been cleaned (solvent/ions removed) and the ligand prepared in MOL2 format.

<div class="grid" markdown>

[:material-download: **protein.pdb** (147 KB)](../assets/examples/md/protein.pdb){ .md-button }

[:material-download: **ligand.mol2** (1.4 KB)](../assets/examples/md/ligand.mol2){ .md-button }

</div>

Or download via command line:

```bash
mkdir prism_basic_tutorial && cd prism_basic_tutorial

wget https://raw.githubusercontent.com/AIB100/PRISM-Tutorial/main/docs/assets/examples/md/protein.pdb
wget https://raw.githubusercontent.com/AIB100/PRISM-Tutorial/main/docs/assets/examples/md/ligand.mol2
```

Verify the files:

```bash
ls -lh protein.pdb ligand.mol2
# protein.pdb  147K  (T4 Lysozyme L99A/M102Q, 164 residues)
# ligand.mol2  1.4K  (Benzene, 12 atoms)
```

## Step 2: Build the System (CLI Method)

### Basic System Building

Use PRISM's command-line interface to build the system:

```bash
# Build with default settings (GAFF force field, amber99sb-ildn, TIP3P water)
prism protein.pdb ligand.mol2 -o t4l_benzene
```

**Expected Output:**
```
🔍 Checking dependencies...
✓ GROMACS found: 2023.3
✓ PDBFixer available
✓ AmberTools found

📋 System Configuration:
  Protein: protein.pdb (164 residues)
  Ligand: ligand.mol2 (12 atoms)
  Output: t4l_benzene
  Ligand FF: GAFF
  Protein FF: amber99sb-ildn
  Water: tip3p

🔧 Building protein-ligand system...
  [1/6] Cleaning protein structure...
  [2/6] Generating ligand parameters (GAFF)...
  [3/6] Building complex topology...
  [4/6] Solvating system...
  [5/6] Adding ions (0.15 M NaCl)...
  [6/6] Generating MDP files...

✅ System built successfully!
📁 Output: ./t4l_benzene/GMX_PROLIG_MD/
```

### Inspect Output Structure

```bash
# Navigate to output
cd t4l_benzene/GMX_PROLIG_MD

# Check generated files
ls -lh

# Expected files:
# - system.gro          # Initial coordinates
# - topol.top           # System topology
# - em/                 # Energy minimization
# - nvt/                # NVT equilibration
# - npt/                # NPT equilibration
# - prod/               # Production run
# - localrun.sh         # Automated run script
```

## Step 3: Run the Simulation

### Option A: Automated Execution

```bash
# Run all stages automatically
bash localrun.sh
```

This script will sequentially run:
1. Energy minimization (EM)
2. NVT equilibration (100 ps)
3. NPT equilibration (100 ps)
4. Production MD (500 ns by default)

### Option B: Step-by-Step Execution

For better control and monitoring:

```bash
# 1. Energy Minimization
cd em
gmx grompp -f em.mdp -c ../system.gro -p ../topol.top -o em.tpr
gmx mdrun -deffnm em -v
cd ..

# Check if EM converged
tail em/em.log
# Look for: "Potential Energy = -XXXXXX.XX"
# Should be large negative number

# 2. NVT Equilibration (constant volume)
cd nvt
gmx grompp -f nvt.mdp -c ../em/em.gro -p ../topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v
cd ..

# 3. NPT Equilibration (constant pressure)
cd npt
gmx grompp -f npt.mdp -c ../nvt/nvt.gro -p ../topol.top -o npt.tpr
gmx mdrun -deffnm npt -v
cd ..

# 4. Production MD
cd prod
gmx grompp -f md.mdp -c ../npt/npt.gro -p ../topol.top -o md.tpr
gmx mdrun -deffnm md -v
cd ..
```

### With GPU Acceleration

```bash
# Single GPU
gmx mdrun -deffnm md -v -nb gpu -pme gpu -bonded gpu -gpu_id 0

# Multiple GPUs
gmx mdrun -deffnm md -v -nb gpu -pme gpu -ntmpi 2 -ntomp 12 -gpu_id 01
```

### Monitoring Progress

```bash
# Check current simulation time
tail -f prod/md.log | grep "Time:"

# Estimate remaining time
gmx check -f prod/md.xtc
```

**Expected Times:**
- CPU only: ~24-48 hours for 500 ns
- Single GPU: ~4-8 hours for 500 ns
- Multi-GPU: ~2-4 hours for 500 ns

## Step 4: Analyze the Trajectory

### Using PRISM's Built-in Analysis

```bash
# Return to project root
cd ../..  # Back to t4l_benzene/

# Run comprehensive analysis
python3 << EOF
import prism as pm

# Analyze trajectory
analysis = pm.analyze_trajectory(
    topology="GMX_PROLIG_MD/system.gro",
    trajectory="GMX_PROLIG_MD/prod/md.xtc",
    ligand_resname="LIG",
    output_dir="analysis_results"
)

print("Analysis complete! Results in analysis_results/")
EOF
```

**Generated Files:**
- `analysis_results/contact_analysis.html` - Interactive contact visualization
- `analysis_results/rmsd_plot.png` - RMSD over time
- `analysis_results/rmsf_plot.png` - Per-residue flexibility
- `analysis_results/contact_heatmap.png` - Protein-ligand contacts
- `analysis_results/distance_analysis.csv` - Detailed distance data

### View Interactive Visualization

```bash
# Open in browser
firefox analysis_results/contact_analysis.html
# or
google-chrome analysis_results/contact_analysis.html
```

### Manual Analysis with GROMACS

Calculate RMSD of protein backbone:

```bash
cd GMX_PROLIG_MD

# Create index groups
gmx make_ndx -f system.gro -o index.ndx << EOF
1 | 13
q
EOF

# Calculate RMSD
echo "Backbone Backbone" | gmx rms -s prod/md.tpr -f prod/md.xtc -o rmsd.xvg -tu ns

# Plot with xmgrace (if installed)
xmgrace rmsd.xvg
```

Calculate protein-ligand distance:

```bash
# Minimum distance between protein and ligand
echo "Protein LIG" | gmx mindist -s prod/md.tpr -f prod/md.xtc -od mindist.xvg -tu ns

# Plot
xmgrace mindist.xvg
```

Calculate hydrogen bonds:

```bash
# Protein-ligand H-bonds
echo "Protein LIG" | gmx hbond -s prod/md.tpr -f prod/md.xtc -num hbond_num.xvg -tu ns

# Plot number of H-bonds over time
xmgrace hbond_num.xvg
```

## Step 5: Visualize the Trajectory

### With PyMOL

```bash
# Load in PyMOL
pymol GMX_PROLIG_MD/system.gro GMX_PROLIG_MD/prod/md.xtc
```

PyMOL commands:
```python
# Align all frames to first frame
intra_fit name CA

# Show protein as cartoon
hide everything
show cartoon, polymer
color slate, polymer

# Show ligand as sticks
show sticks, resn LIG
color green, resn LIG
util.cbag resn LIG

# Play trajectory
mplay
```

### With VMD

```bash
vmd GMX_PROLIG_MD/system.gro GMX_PROLIG_MD/prod/md.xtc
```

VMD console commands:
```tcl
# Align trajectory
set sel [atomselect top "protein and name CA"]
$sel frame 0
set ref [atomselect top "protein and name CA" frame 0]

for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
    $sel frame $i
    set trans [measure fit $sel $ref]
    [atomselect top all frame $i] move $trans
}

# Representations
mol delrep 0 top
mol representation NewCartoon
mol selection protein
mol addrep top

mol representation Licorice
mol selection resname LIG
mol addrep top
```

## Step 6: Interpret Results

### Expected Observations

For a well-behaved simulation, you should see:

1. **RMSD Stabilization**
   - Protein RMSD should plateau after ~10-50 ns
   - Typical values: 1-3 Å for backbone
   - Ligand RMSD may fluctuate more (normal)

2. **Stable Protein-Ligand Contacts**
   - Consistent hydrophobic contacts
   - Occasional hydrogen bonds (if capable)
   - Ligand stays in binding pocket

3. **Reasonable Energies**
   - Potential energy: Large negative value (-4×10⁵ to -6×10⁵ kJ/mol for this system)
   - Temperature: Should average ~310 K
   - Pressure: Should average ~1 bar

### Quality Checks

```bash
# Check temperature
echo "Temperature" | gmx energy -f prod/md.edr -o temperature.xvg

# Check pressure
echo "Pressure" | gmx energy -f prod/md.edr -o pressure.xvg

# Check potential energy
echo "Potential" | gmx energy -f prod/md.edr -o potential.xvg

# Check total energy
echo "Total-Energy" | gmx energy -f prod/md.edr -o total_energy.xvg
```

## Troubleshooting

### Simulation Crashes

**"LINCS WARNING"**
- System unstable, likely due to clashes
- Solution: Re-run energy minimization with stricter parameters
- Or: Reduce time step to 1 fs

**"Segmentation fault"**
- Usually GPU-related
- Try: Run on CPU only first
- Check: GPU drivers and CUDA version

### Analysis Errors

**"No frames found"**
- Simulation may not have run
- Check: `ls -lh prod/md.xtc` should show file size

**"Residue LIG not found"**
- Ligand residue name mismatch
- Check: `grep "resname" system.gro` to find actual residue name

## Next Steps

Congratulations! You've completed the basic PRISM tutorial. Here's where to go next:

1. **Try Different Force Fields**
   - Rebuild with OpenFF: `prism protein.pdb ligand.mol2 -o openff_system --ligand-forcefield openff`
   - Compare results with GAFF version
   - See [Force Field Tutorial](force-field-tutorial.md)

2. **Calculate Binding Energy**
   - Run PMF calculations on your equilibrated system
   - See [PMF Tutorial](pmf-tutorial.md)

3. **Process Multiple Ligands**
   - Automate building for virtual screening
   - See [Batch Tutorial](batch-tutorial.md)

4. **Advanced Analysis**
   - Free energy decomposition
   - Contact network analysis
   - See [Analysis Tools Guide](../user-guide/analysis-tools.md)

## Summary

In this tutorial, you learned:

✅ How to prepare protein and ligand inputs
✅ Building MD systems with PRISM CLI
✅ Running equilibration and production simulations
✅ Analyzing protein-ligand interactions
✅ Visualizing MD trajectories
✅ Interpreting simulation quality

## Additional Resources

- [PRISM User Guide](../user-guide/index.md)
- [Force Fields Guide](../user-guide/force-fields.md)
- [Running Simulations](../user-guide/running-simulations.md)
- [Analysis Tools](../user-guide/analysis-tools.md)

<div class="whats-next" markdown>

## What's Next

- [Compare force fields in the Force Field Tutorial](force-field-tutorial.md)
- [Calculate binding energy with the PMF Tutorial](pmf-tutorial.md)
- [Automate builds with the Batch Processing Tutorial](batch-tutorial.md)
- [Explore the full User Guide](../user-guide/index.md)

</div>
