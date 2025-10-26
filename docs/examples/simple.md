# Simple Protein-Ligand Example

A complete, minimal example demonstrating the PRISM workflow from protein and ligand files to a production-ready MD system.

## Overview

This example shows how to:

- Build a protein-ligand system with default settings
- Generate all necessary GROMACS files
- Run molecular dynamics simulations
- Perform basic trajectory analysis

**Complexity**: Beginner  
**Time**: 30 minutes + simulation time  
**System**: T4 Lysozyme + Benzene

---

## Prerequisites

### Required Files

Download the example files:

```bash
# Create working directory
mkdir simple_example
cd simple_example

# Download protein (T4 Lysozyme)
wget https://files.rcsb.org/download/2LZM.pdb -O protein.pdb

# Create benzene ligand (simplified MOL2 format)
cat > ligand.mol2 << 'LIGAND'
@<TRIPOS>MOLECULE
BNZ
 12 12 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 C1          0.0000    1.3970    0.0000 C.ar    1  BNZ1       -0.1150
      2 C2          1.2098    0.6985    0.0000 C.ar    1  BNZ1       -0.1150
      3 C3          1.2098   -0.6985    0.0000 C.ar    1  BNZ1       -0.1150
      4 C4          0.0000   -1.3970    0.0000 C.ar    1  BNZ1       -0.1150
      5 C5         -1.2098   -0.6985    0.0000 C.ar    1  BNZ1       -0.1150
      6 C6         -1.2098    0.6985    0.0000 C.ar    1  BNZ1       -0.1150
      7 H1          0.0000    2.4810    0.0000 H       1  BNZ1        0.1150
      8 H2          2.1488    1.2405    0.0000 H       1  BNZ1        0.1150
      9 H3          2.1488   -1.2405    0.0000 H       1  BNZ1        0.1150
     10 H4          0.0000   -2.4810    0.0000 H       1  BNZ1        0.1150
     11 H5         -2.1488   -1.2405    0.0000 H       1  BNZ1        0.1150
     12 H6         -2.1488    1.2405    0.0000 H       1  BNZ1        0.1150
@<TRIPOS>BOND
     1     1     2   ar
     2     2     3   ar
     3     3     4   ar
     4     4     5   ar
     5     5     6   ar
     6     6     1   ar
     7     1     7    1
     8     2     8    1
     9     3     9    1
    10     4    10    1
    11     5    11    1
    12     6    12    1
@<TRIPOS>SUBSTRUCTURE
     1 BNZ1        1 TEMP              0 ****  ****    0 ROOT
LIGAND
```

Your directory should now contain:
```
simple_example/
├── protein.pdb
└── ligand.mol2
```

---

## Complete Script

Create `run_simple_example.py`:

```python
#!/usr/bin/env python3
"""
Simple PRISM Example
Build a basic protein-ligand MD system
"""

import prism
import os

def main():
    print("="*60)
    print("PRISM Simple Example: T4 Lysozyme + Benzene")
    print("="*60)
    
    # Input files
    protein = "protein.pdb"
    ligand = "ligand.mol2"
    output_dir = "t4_lysozyme_benzene"
    
    # Check input files exist
    if not os.path.exists(protein):
        print(f"ERROR: {protein} not found!")
        return
    if not os.path.exists(ligand):
        print(f"ERROR: {ligand} not found!")
        return
    
    print(f"\nInput files:")
    print(f"  Protein: {protein}")
    print(f"  Ligand:  {ligand}")
    print(f"  Output:  {output_dir}")
    
    # Create system
    print("\nStep 1: Creating PRISM system...")
    system = prism.system(
        protein,
        ligand,
        output_dir=output_dir
    )
    
    # Build system (this does everything!)
    print("\nStep 2: Building system...")
    print("This will:")
    print("  - Generate ligand force field parameters (GAFF)")
    print("  - Clean and prepare protein structure")
    print("  - Build GROMACS topology")
    print("  - Solvate system in water box")
    print("  - Add ions for neutralization and 0.15 M salt")
    print("  - Generate MDP files for EM, NVT, NPT, Production")
    print("  - Create ready-to-run simulation script")
    print()
    
    system.build()
    
    print("\n" + "="*60)
    print("SUCCESS! System built successfully")
    print("="*60)
    
    # Display output structure
    print("\nGenerated files:")
    print(f"{output_dir}/")
    print("├── LIG.amb2gmx/          # Ligand force field")
    print("│   ├── LIG.gro          # Ligand coordinates")
    print("│   ├── LIG.itp          # Ligand topology")
    print("│   └── ...")
    print("├── GMX_PROLIG_MD/        # Main simulation directory")
    print("│   ├── solv_ions.gro    # Complete solvated system")
    print("│   ├── topol.top        # System topology")
    print("│   └── localrun.sh      # Run script")
    print("├── mdps/                 # MD parameter files")
    print("│   ├── em.mdp")
    print("│   ├── nvt.mdp")
    print("│   ├── npt.mdp")
    print("│   └── md.mdp")
    print("└── prism_config.yaml     # Configuration used")
    
    print("\nNext steps:")
    print(f"  1. cd {output_dir}/GMX_PROLIG_MD")
    print("  2. bash localrun.sh")
    print("  3. Wait for simulation to complete (~2-4 hours with GPU)")
    print()

if __name__ == "__main__":
    main()
```

---

## Running the Example

### Step 1: Build the System

```bash
python run_simple_example.py
```

Expected output:
```
============================================================
PRISM Simple Example: T4 Lysozyme + Benzene
============================================================

Input files:
  Protein: protein.pdb
  Ligand:  ligand.mol2
  Output:  t4_lysozyme_benzene

Step 1: Creating PRISM system...

Step 2: Building system...
This will:
  - Generate ligand force field parameters (GAFF)
  - Clean and prepare protein structure
  - Build GROMACS topology
  - Solvate system in water box
  - Add ions for neutralization and 0.15 M salt
  - Generate MDP files for EM, NVT, NPT, Production
  - Create ready-to-run simulation script

[PRISM] Generating ligand force field...
[PRISM] Cleaning protein structure...
[PRISM] Building GROMACS model...
[PRISM] Adding solvent and ions...
[PRISM] Generating MDP files...
[PRISM] Creating run script...

============================================================
SUCCESS! System built successfully
============================================================
```

### Step 2: Run MD Simulation

```bash
cd t4_lysozyme_benzene/GMX_PROLIG_MD
bash localrun.sh
```

The script will automatically run:
1. Energy Minimization (~5 minutes)
2. NVT Equilibration (~10 minutes)
3. NPT Equilibration (~10 minutes)
4. Production MD (~2-4 hours for 500 ns)

### Step 3: Monitor Progress

While simulation is running:

```bash
# Check current stage
ls -lh prod/

# Monitor energy
tail -f prod/md.log

# Check if finished
ls prod/md.xtc
```

---

## Analyzing Results

After simulation completes:

```python
#!/usr/bin/env python3
"""Analyze the simulation results"""

import prism

# Analyze trajectory
analyzer = prism.analyze_trajectory(
    topology="t4_lysozyme_benzene/GMX_PROLIG_MD/prod/md.tpr",
    trajectory="t4_lysozyme_benzene/GMX_PROLIG_MD/prod/md.xtc",
    ligand_resname="LIG",
    output_dir="analysis_results"
)

print("Analysis complete! Check analysis_results/ directory")
```

This generates:
- RMSD plots (protein backbone and ligand)
- RMSF per-residue fluctuations
- Contact frequency maps
- Hydrogen bond analysis
- Interactive HTML visualization

---

## Command-Line Alternative

You can also use the command-line interface:

```bash
# Build system
prism protein.pdb ligand.mol2 -o t4_lysozyme_benzene

# Run simulation
cd t4_lysozyme_benzene/GMX_PROLIG_MD
bash localrun.sh
```

---

## Expected Results

### System Statistics

After building, you should have:
- Total atoms: ~40,000-60,000
- Protein atoms: ~2,500
- Ligand atoms: 12
- Water molecules: ~12,000-18,000
- Na+ and Cl- ions: ~30-50

### Simulation Metrics

Expected values for T4 Lysozyme-Benzene:
- Protein backbone RMSD: 1-2 Å
- Ligand RMSD: 1-3 Å (if bound)
- Protein-ligand contacts: 5-15
- Hydrogen bonds: 0-2 (benzene has no donors/acceptors)

---

## Customizing the Example

### Change Force Field

```python
system = prism.system(
    protein, ligand,
    output_dir=output_dir,
    ligand_forcefield="openff",  # Use OpenFF instead of GAFF
    forcefield="amber14sb"        # Use AMBER14SB for protein
)
```

### Shorter Simulation

```python
system = prism.system(
    protein, ligand,
    output_dir=output_dir,
    production_ns=100  # 100 ns instead of 500 ns
)
```

### Different Water Model

```python
system = prism.system(
    protein, ligand,
    output_dir=output_dir,
    water_model="tip4p"  # TIP4P water
)
```

---

## Troubleshooting

### Issue: "GROMACS command not found"

**Solution**: Source GROMACS
```bash
source /path/to/gromacs/bin/GMXRC
gmx --version
```

### Issue: "Cannot generate GAFF parameters"

**Solution**: Install AmberTools
```bash
conda install -c conda-forge ambertools
pip install acpype
```

### Issue: "Ligand parameterization failed"

**Solution**: Try a different force field
```python
ligand_forcefield="openff"  # More robust
```

---

## Next Examples

- [Multiple Force Fields](multi-forcefield.md) - Compare different force fields
- [Custom Configuration](custom-config.md) - Advanced customization
- [Batch Processing](../tutorials/batch-tutorial.md) - Process multiple ligands

---

## Complete Working Example

Download the complete example:

```bash
git clone https://github.com/AIB001/PRISM-examples.git
cd PRISM-examples/simple-protein-ligand
python run_example.py
```

---

## Questions?

- See [Troubleshooting Guide](../user-guide/troubleshooting.md)
- Check [API Documentation](../api/index.md)
- Ask in [GitHub Discussions](https://github.com/AIB001/PRISM/discussions)
