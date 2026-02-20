# Quick Start Guide

Get started with PRISM in 5 minutes! This guide demonstrates the complete workflow from preparing input files to running MD simulations.

!!! example "Quick Start"
    ```bash
    prism protein.pdb ligand.mol2 -o my_simulation
    cd my_simulation/GMX_PROLIG_MD && bash localrun.sh
    ```

## Prerequisites Check

Before starting, ensure you have:
- ✅ PRISM installed (`pip install -e .`)
- ✅ GROMACS 2024.3+ installed and accessible
- ✅ Force field dependencies installed (see [Installation](installation.md))
- ✅ Your protein PDB file (with or without hydrogens)
- ✅ Your ligand MOL2/SDF file (with correct protonation state)

## Method 1: Command-Line Interface (CLI)

### Step 1: Build the System

The simplest way - build with one command:

```bash
# Default: GAFF force field, amber99sb protein, tip3p water
prism protein.pdb ligand.mol2 -o my_simulation

# The system will be built in ./my_simulation directory
```

What PRISM does automatically:
1. ✓ Cleans protein structure (removes HETATM, handles metals)
2. ✓ Generates ligand force field parameters (GAFF by default)
3. ✓ Builds GROMACS topology files
4. ✓ Solvates system in water box
5. ✓ Neutralizes and adds physiological salt (0.15 M)
6. ✓ Generates MDP files for all simulation stages
7. ✓ Creates ready-to-run simulation script

### Step 2: Run MD Simulation

After building, run the simulation:

```bash
cd my_simulation/GMX_PROLIG_MD
bash localrun.sh
```

The script will automatically run:
1. Energy Minimization (EM)
2. NVT Equilibration (temperature)
3. NPT Equilibration (pressure)
4. Production MD (data collection)

!!! success "That's it!"
    Your system is now running molecular dynamics. Results will be in `prod/` directory.

---

## Method 2: Python API

For more control, use the Python interface:

```python
import prism

# Create system
system = prism.system("protein.pdb", "ligand.mol2", output_dir="my_sim")

# Build the system
system.build()

# Check what was generated
system.info()

# Get file paths
files = system.get_output_files()
print(files['system_gro'])  # Path to solvated system
print(files['system_top'])  # Path to topology
```

---

## Choosing Force Fields

### Using Different Ligand Force Fields

PRISM supports 8+ ligand force fields:

```bash
# GAFF (default, most common)
prism protein.pdb ligand.mol2 -o output

# GAFF2 (improved version)
prism protein.pdb ligand.mol2 -o output --ligand-forcefield gaff2

# OpenFF (modern, data-driven)
prism protein.pdb ligand.sdf -o output --ligand-forcefield openff

# OPLS-AA (via LigParGen server)
prism protein.pdb ligand.mol2 -o output --ligand-forcefield opls

# CGenFF (CHARMM, requires web download)
prism protein.pdb ligand.mol2 -o output --ligand-forcefield cgenff \
  --forcefield-path /path/to/downloaded_cgenff

# MMFF/MATCH/Hybrid (via SwissParam)
prism protein.pdb ligand.mol2 -o output --ligand-forcefield mmff
```

### Using Different Protein Force Fields

```bash
# AMBER14SB (recommended for proteins)
prism protein.pdb ligand.mol2 -o output --forcefield amber14sb

# AMBER99SB-ILDN (improved side chains)
prism protein.pdb ligand.mol2 -o output --forcefield amber99sb-ildn

# CHARMM36 (good for membranes)
prism protein.pdb ligand.mol2 -o output --forcefield charmm36 --water tips3p
```

---

## Understanding the Output

After PRISM completes, you'll find this directory structure:

```
my_simulation/
├── LIG.amb2gmx/              # Ligand force field files (GAFF/GAFF2)
│   ├── LIG.gro              # Ligand coordinates
│   ├── LIG.itp              # Ligand topology
│   ├── atomtypes_LIG.itp    # Atom type definitions
│   └── posre_LIG.itp        # Position restraints
│
├── LIG.openff2gmx/           # Ligand files (if using OpenFF)
│   └── (same structure)
│
├── GMX_PROLIG_MD/            # **Main simulation directory**
│   ├── solv_ions.gro        # Complete solvated system
│   ├── topol.top            # System topology
│   ├── localrun.sh          # Ready-to-run simulation script
│   └── (em/, nvt/, npt/, prod/ created during simulation)
│
├── mdps/                     # MD parameter files
│   ├── em.mdp               # Energy minimization
│   ├── nvt.mdp              # NVT equilibration (500 ps)
│   ├── npt.mdp              # NPT equilibration (500 ps)
│   └── md.mdp               # Production MD (500 ns default)
│
├── prism_config.yaml         # Configuration used
└── protein_clean.pdb         # Cleaned protein structure
```

### Key Files Explained

| File | Purpose |
|------|---------|
| `solv_ions.gro` | Starting structure for MD (protein + ligand + water + ions) |
| `topol.top` | Complete system topology (defines all interactions) |
| `localrun.sh` | Automated script to run complete MD workflow |
| `*.mdp` files | Parameters for each simulation stage |
| `LIG.itp` | Ligand-specific parameters (charges, bonds, angles, dihedrals) |

---

## Advanced Options

### Customizing Simulation Parameters

```bash
# Longer production run (1 microsecond)
prism protein.pdb ligand.mol2 -o output --production-ns 1000

# Different temperature (300 K instead of 310 K)
prism protein.pdb ligand.mol2 -o output --temperature 300

# Larger water box (2.0 nm instead of 1.5 nm)
prism protein.pdb ligand.mol2 -o output --box-distance 2.0

# Higher salt concentration (physiological is 0.15 M)
prism protein.pdb ligand.mol2 -o output --salt-concentration 0.15

# Charged ligand
prism protein.pdb ligand.mol2 -o output --ligand-charge -1
```

### Using Configuration Files

For complex setups, use a YAML configuration file:

```bash
# Export default configuration template
prism --export-config my_config.yaml

# Edit my_config.yaml with your preferred settings

# Build with custom config
prism protein.pdb ligand.mol2 -o output --config my_config.yaml
```

Example `my_config.yaml`:
```yaml
simulation:
  temperature: 300
  pressure: 1.0
  pH: 7.4
  production_time_ns: 1000

box:
  distance: 2.0
  shape: dodecahedron  # More efficient than cubic

ions:
  concentration: 0.15  # Physiological salt
```

---

## Common Workflows

### Workflow 1: Drug-Protein Complex

```bash
# Download PDB from RCSB
# Prepare ligand with correct protonation (use ChemDraw, Avogadro, etc.)

# Build system with default GAFF
prism 1ABC_protein.pdb drug.mol2 -o drug_complex

# Run MD
cd drug_complex/GMX_PROLIG_MD
bash localrun.sh
```

### Workflow 2: High-Accuracy Binding Study

```bash
# Use OpenFF for better ligand parameters
# Use AMBER14SB for better protein backbone
prism receptor.pdb compound.sdf -o high_accuracy \
  --ligand-forcefield openff \
  --forcefield amber14sb \
  --production-ns 1000 \
  --box-distance 2.0
```

### Workflow 3: Screening Multiple Ligands

```python
import prism
from pathlib import Path

# List of ligand files
ligands = list(Path("ligands/").glob("*.mol2"))

for lig in ligands:
    name = lig.stem
    system = prism.system("protein.pdb", str(lig), output_dir=f"screen/{name}")
    try:
        system.build()
        print(f"✓ Built system for {name}")
    except Exception as e:
        print(f"✗ Failed for {name}: {e}")
```

---

## Analysis with PRISM

PRISM includes comprehensive analysis tools:

### Quick Analysis

```python
import prism

# Analyze trajectory and generate interactive visualization
analyzer = prism.analyze_trajectory(
    topology="system.gro",
    trajectory="md.xtc",
    ligand_resname="LIG"
)

# Results saved to analysis_results/ directory
# Includes:
# - RMSD/RMSF plots
# - Contact frequency maps
# - Hydrogen bond analysis
# - Interactive HTML visualization
```

### Interactive Contact Visualization

```python
import prism

# Generate interactive HTML contact map
prism.visualize_trajectory(
    trajectory="prod/md.xtc",
    topology="prod/md.tpr",
    ligand="ligand.sdf",
    output="contacts.html"
)

# Open contacts.html in web browser for interactive viewing
```

### Basic GROMACS Analysis

```bash
cd my_simulation/GMX_PROLIG_MD/prod

# RMSD (backbone stability)
echo "Backbone Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg

# RMSF (per-residue fluctuations)
echo "C-alpha" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg

# Protein-ligand distance
echo "Protein LIG" | gmx distance -s md.tpr -f md.xtc -oall distances.xvg

# Hydrogen bonds
echo "Protein LIG" | gmx hbond -s md.tpr -f md.xtc -num hbonds.xvg
```

---

## Command-Line Options Reference

| Option | Description | Default | Example |
|--------|-------------|---------|---------|
| `-o, --output` | Output directory | `prism_output` | `-o my_system` |
| `--ligand-forcefield` | Ligand force field | `gaff` | `--ligand-forcefield openff` |
| `--forcefield` | Protein force field | `amber99sb` | `--forcefield amber14sb` |
| `--water` | Water model | `tip3p` | `--water tip4p` |
| `--box-distance` | Box padding (nm) | `1.5` | `--box-distance 2.0` |
| `--temperature` | Temperature (K) | `310` | `--temperature 300` |
| `--production-ns` | Production time (ns) | `500` | `--production-ns 1000` |
| `--salt-concentration` | Salt conc. (M) | `0.15` | `--salt-concentration 0.10` |
| `--ligand-charge` | Ligand net charge | `0` | `--ligand-charge -1` |
| `--config` | Config YAML file | None | `--config custom.yaml` |
| `--overwrite` | Overwrite existing | `False` | `--overwrite` |

---

## Troubleshooting

### Common Build Issues

!!! failure "GROMACS not found"
    ```
    Error: GROMACS command not found
    ```
    **Solution**: Ensure GROMACS is installed and sourced:
    ```bash
    source /path/to/gromacs/bin/GMXRC
    gmx --version
    ```

!!! failure "Force field not available"
    ```
    Error: Cannot generate GAFF parameters
    ```
    **Solution**: Install AmberTools and ACPYPE:
    ```bash
    conda install -c conda-forge ambertools
    pip install acpype
    ```

!!! failure "Ligand parameterization failed"
    ```
    Error: Cannot assign atom types
    ```
    **Solutions**:
    - Check ligand file format (MOL2 should have correct atom types)
    - Try a different force field: `--ligand-forcefield openff`
    - Verify ligand has proper 3D coordinates and hydrogens

### Common Simulation Issues

!!! failure "System has non-zero charge"
    **Solution**: Check ligand charge and specify if needed:
    ```bash
    prism protein.pdb ligand.mol2 -o output --ligand-charge -1
    ```

!!! failure "Energy minimization fails"
    **Solution**: Increase box size to reduce clashes:
    ```bash
    prism protein.pdb ligand.mol2 -o output --box-distance 2.0
    ```

!!! failure "GPU not detected"
    **Solution**:
    - Check: `nvidia-smi` and `gmx mdrun -version`
    - Edit `localrun.sh` to remove `-gpu_id 0` flags for CPU-only

---

<div class="whats-next" markdown>

## What's Next

- [Explore force field options in the Force Fields Guide](../user-guide/force-fields.md)
- [Learn about analysis tools for your trajectories](../user-guide/analysis-tools.md)
- [Set up PMF calculations for binding free energy](../user-guide/pmf-calculations.md)
- [Read the full User Guide for detailed documentation](../user-guide/index.md)

</div>