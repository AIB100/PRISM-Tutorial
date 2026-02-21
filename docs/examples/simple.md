# Simple Protein-Ligand Example

!!! example "Quick Start"
    ```bash
    prism protein.pdb ligand.mol2 -o t4l_benzene
    cd t4l_benzene/GMX_PROLIG_MD && bash localrun.sh
    ```

A complete, minimal example demonstrating the PRISM workflow from input files to a production-ready MD system.

## Overview

This example shows how to:

- Build a protein-ligand system with default settings
- Generate all necessary GROMACS files
- Run molecular dynamics simulations

**Complexity**: Beginner
**Time**: ~10 minutes (system building) + simulation time
**System**: T4 Lysozyme L99A/M102Q + Benzene (PDB: [5JWT](https://www.rcsb.org/structure/5JWT))

---

## Download Example Files

Pre-processed input files are provided. The protein structure has been cleaned (solvent/ions removed, hydrogens added) and the ligand has been prepared in MOL2 format with correct atom types.

<div class="grid" markdown>

[:material-download: **protein.pdb** (112 KB)](../assets/examples/md/protein.pdb){ .md-button }

[:material-download: **ligand.mol2** (1 KB)](../assets/examples/md/ligand.mol2){ .md-button }

</div>

Or download via command line:

```bash
mkdir simple_example && cd simple_example

# Download from the PRISM Tutorial repository
wget https://raw.githubusercontent.com/AIB100/PRISM-Tutorial/main/docs/assets/examples/md/protein.pdb
wget https://raw.githubusercontent.com/AIB100/PRISM-Tutorial/main/docs/assets/examples/md/ligand.mol2
```

Your directory should contain:

```
simple_example/
├── protein.pdb    # T4 Lysozyme L99A/M102Q (164 residues, 1309 atoms)
└── ligand.mol2    # Benzene (12 atoms: 6C + 6H)
```

---

## Step 1: Build the System

### Command Line

```bash
prism protein.pdb ligand.mol2 -o t4l_benzene
```

### Python API

```python
import prism as pm

system = pm.system("protein.pdb", "ligand.mol2", output_dir="t4l_benzene")
system.build()
```

PRISM will automatically:

1. Generate ligand force field parameters (GAFF by default)
2. Prepare the protein topology (AMBER99SB-ILDN by default)
3. Combine protein and ligand into a complex
4. Solvate in a TIP3P water box
5. Add Na+/Cl- ions for charge neutralization (0.15 M)
6. Generate MDP files for all simulation stages
7. Create a ready-to-run simulation script

### Output Structure

```
t4l_benzene/
├── LIG.amb2gmx/              # Ligand force field files
│   ├── LIG.gro               # Ligand coordinates
│   ├── LIG.itp               # Ligand topology
│   └── ...
├── GMX_PROLIG_MD/            # Simulation directory
│   ├── solv_ions.gro         # Solvated system
│   ├── topol.top             # System topology
│   └── localrun.sh           # Run script
└── mdps/                     # MD parameter files
    ├── em.mdp                # Energy minimization
    ├── nvt.mdp               # NVT equilibration
    ├── npt.mdp               # NPT equilibration
    └── md.mdp                # Production MD
```

---

## Step 2: Run the Simulation

```bash
cd t4l_benzene/GMX_PROLIG_MD
bash localrun.sh
```

The script runs four stages sequentially:

| Stage | Purpose | Typical Time (GPU) |
| --- | --- | --- |
| Energy Minimization | Remove bad contacts | ~2 minutes |
| NVT Equilibration | Heat to 300 K | ~10 minutes |
| NPT Equilibration | Equilibrate pressure | ~10 minutes |
| Production MD | Data collection | 2-4 hours (for 100 ns) |

### Monitor Progress

```bash
# Check which stage is running
ls -lh em/ nvt/ npt/ prod/

# Watch the production MD log
tail -f prod/md.log
```

---

## Step 3: Check Results

After simulation completes, verify the output:

```bash
# Check production trajectory exists
ls -lh prod/md.xtc

# Quick energy check
gmx energy -f prod/md.edr -o energy.xvg
```

### Expected System Statistics

| Property | Expected Value |
| --- | --- |
| Total atoms | ~25,000-35,000 |
| Protein atoms | 1,309 |
| Ligand atoms | 12 |
| Water molecules | ~8,000-11,000 |
| Na+/Cl- ions | ~20-40 |

---

## Customization

### Change Force Field

```bash
# Use OpenFF for ligand, AMBER14SB for protein
prism protein.pdb ligand.mol2 -o t4l_benzene \
    --ligand-forcefield openff --forcefield amber14sb
```

### Adjust Simulation Length

Edit `mdps/md.mdp` before running, or use configuration:

```bash
prism protein.pdb ligand.mol2 -o t4l_benzene --production-time 500
```

### Change Water Model

```bash
prism protein.pdb ligand.mol2 -o t4l_benzene --water tip4p
```

---

## Troubleshooting

**"GROMACS command not found"**: Source the GROMACS environment first:
```bash
source /path/to/gromacs/bin/GMXRC
```

**"Cannot generate GAFF parameters"**: Install AmberTools and ACPYPE:
```bash
conda install -c conda-forge ambertools
pip install acpype
```

**"Ligand parameterization failed"**: Try OpenFF as an alternative:
```bash
prism protein.pdb ligand.mol2 -o t4l_benzene --ligand-forcefield openff
```

---

<div class="whats-next" markdown>

## What's Next

- [Compare force fields in the Multi-Force Field Example](multi-forcefield.md)
- [Customize parameters in the Custom Config Example](custom-config.md)
- [Submit to a cluster with Slurm Templates](../tutorials/cluster-submission.md)

</div>
