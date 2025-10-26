# Force Field Selection Tutorial

## Overview

Learn how to select and compare different force fields for your protein-ligand system using PRISM's support for 8+ ligand parameterization methods.

**Time:** 45-90 minutes
**Level:** Intermediate

## Prerequisites

- Completed [Basic Tutorial](basic-tutorial.md)
- PRISM with force field dependencies installed
- Understanding of force field concepts

## Force Fields Covered

This tutorial demonstrates:

- **GAFF** - General AMBER Force Field
- **GAFF2** - Improved GAFF
- **OpenFF** - Open Force Field
- **OPLS-AA** - Using LigParGen

## Step 1: Prepare the System

```bash
mkdir ff_comparison
cd ff_comparison

# Use same inputs as basic tutorial
wget https://files.rcsb.org/download/3DMX.pdb -O protein.pdb
wget https://github.com/AIB001/PRISM-tutorial-data/raw/main/basic/benzene.mol2
```

## Step 2: Build with Different Force Fields

### GAFF (Default)

```bash
prism protein.pdb benzene.mol2 -o gaff_system --ligand-ff gaff
```

### GAFF2

```bash
prism protein.pdb benzene.mol2 -o gaff2_system --ligand-ff gaff2
```

### OpenFF

```bash
prism protein.pdb benzene.mol2 -o openff_system --ligand-ff openff
```

### OPLS-AA

```bash
prism protein.pdb benzene.mol2 -o opls_system --ligand-ff opls
```

## Step 3: Run Simulations

```bash
# Run all systems (can be parallelized)
for system in gaff_system gaff2_system openff_system opls_system; do
    cd $system/GMX_PROLIG_MD
    bash localrun.sh &
    cd ../..
done
```

## Step 4: Compare Results

### Binding Stability

```python
import prism as pm
import matplotlib.pyplot as plt
import numpy as np

force_fields = ['gaff', 'gaff2', 'openff', 'opls']
results = {}

for ff in force_fields:
    analysis = pm.analyze_trajectory(
        f"{ff}_system/GMX_PROLIG_MD/system.gro",
        f"{ff}_system/GMX_PROLIG_MD/prod/md.xtc",
        ligand_resname="BEN",
        output_dir=f"{ff}_analysis"
    )
    results[ff] = analysis

# Compare RMSD, contacts, etc.
print("Force Field Comparison:")
for ff in force_fields:
    print(f"{ff.upper()}: {results[ff]['summary']}")
```

## Step 5: Interpret Differences

Compare force fields based on:

1. **Charge Assignment**
   - Check partial charges in topology files
   - Different methods yield different charges

2. **Binding Stability**
   - Which FF keeps ligand most stable?
   - RMSD fluctuations

3. **Interaction Patterns**
   - Hydrogen bonding differences
   - Hydrophobic contacts

## Recommendations

- **GAFF/GAFF2**: Good default, widely tested
- **OpenFF**: More accurate for drug-like molecules
- **OPLS-AA**: Excellent for aromatic systems
- **CGenFF**: Best for CHARMM protein FF

## Next Steps

- [PMF Tutorial](pmf-tutorial.md) - Calculate binding energies
- [User Guide: Force Fields](../user-guide/force-fields.md)

