# PMF Calculation Example

!!! example "Quick Start"
    ```bash
    prism protein.pdb ligand.mol2 -ff amber14sb -lff gaff2 \
        --gaussian hf --isopt false -o CDK2-CVT313 --pmf
    ```

A complete example of setting up a PMF (Potential of Mean Force) calculation to estimate the binding free energy of a protein-ligand complex.

## Overview

This example demonstrates:

- Building a PMF-ready system with optimized pulling direction
- Using Gaussian HF charges for accurate ligand electrostatics
- Running steered MD, umbrella sampling, and WHAM analysis
- Extracting the binding free energy from the PMF profile

**Complexity**: Intermediate
**Time**: ~15 minutes (system building) + simulation time
**System**: CDK2 + CVT-313 inhibitor (PDB: [6INL](https://www.rcsb.org/structure/6INL))

---

## About the System

**Cyclin-dependent kinase 2 (CDK2)** is a key cell cycle regulator and an important drug target in oncology. **CVT-313** is a purine-based inhibitor that binds to the ATP-binding pocket of CDK2 with nanomolar affinity.

The original crystal structure (PDB: 6INL) has missing residues. For this example, the protein was completed using **Swiss-Model** homology modeling based on the crystal structure, ensuring a continuous backbone suitable for MD simulations.

---

## Download Example Files

<div class="grid" markdown>

[:material-download: **protein.pdb** (197 KB)](../assets/examples/pmf/protein.pdb){ .md-button }

[:material-download: **ligand.mol2** (3.6 KB)](../assets/examples/pmf/ligand.mol2){ .md-button }

</div>

Or download via command line:

```bash
mkdir pmf_example && cd pmf_example

wget https://raw.githubusercontent.com/AIB100/PRISM-Tutorial/main/docs/assets/examples/pmf/protein.pdb
wget https://raw.githubusercontent.com/AIB100/PRISM-Tutorial/main/docs/assets/examples/pmf/ligand.mol2
```

| File | Description |
| --- | --- |
| `protein.pdb` | CDK2, Swiss-Model homology model from 6INL (299 residues, 2398 atoms) |
| `ligand.mol2` | CVT-313 inhibitor (57 atoms, purine scaffold with phenyl substituent) |

---

## Step 1: Build the PMF System

```bash
prism protein.pdb ligand.mol2 \
    -ff amber14sb \
    -lff gaff2 \
    --gaussian hf \
    --isopt false \
    -o CDK2-CVT313 \
    --pmf
```

### Flag Explanation

| Flag | Purpose |
| --- | --- |
| `-ff amber14sb` | Use AMBER14SB force field for the protein |
| `-lff gaff2` | Use GAFF2 force field for the ligand |
| `--gaussian hf` | Compute RESP charges at HF/6-31G* level using Gaussian |
| `--isopt false` | Skip geometry optimization (use crystal geometry directly) |
| `-o CDK2-CVT313` | Output directory name |
| `--pmf` | Enable PMF mode: align system, extend box, generate SMD/umbrella files |

!!! note "Why Gaussian RESP charges?"
    CVT-313 is a multi-ring heterocyclic compound. Standard AM1-BCC charges may not adequately capture the electrostatic profile of such systems. HF/6-31G* RESP charges provide more accurate partial charges, especially for nitrogen-containing heterocycles.

!!! note "Why skip geometry optimization?"
    The ligand coordinates come from a crystal structure (experimental geometry). For molecules with well-determined X-ray geometry, optimization is unnecessary and skipping it saves computational time.

### What PRISM Does Automatically

1. **Gaussian charge calculation**: Submit HF/6-31G* ESP calculation, fit RESP charges
2. **GAFF2 parameterization**: Generate ligand topology with RESP charges
3. **System assembly**: Combine protein + ligand, solvate, add ions
4. **Pulling direction optimization**: Metropolis-Hastings simulated annealing on the unit sphere to find the least-obstructed pulling path
5. **System alignment**: Rotate the complex so the optimal pulling direction aligns with +Z
6. **Box extension**: Extend the simulation box along Z to accommodate ligand unbinding
7. **MDP generation**: Create parameter files for SMD and umbrella sampling

### Output Structure

```
CDK2-CVT313/
├── LIG.amb2gmx/                # Ligand GAFF2 parameters (with RESP charges)
├── GMX_PROLIG_MD/              # Simulation directory
│   ├── solv_ions.gro           # Aligned, solvated system
│   ├── topol.top               # System topology
│   ├── run_smd.sh              # Steered MD script
│   ├── run_umbrella.sh         # Umbrella sampling script
│   ├── run_wham.sh             # WHAM analysis script
│   └── localrun.sh             # Sequential all-stages script
├── mdps/
│   ├── em.mdp                  # Energy minimization
│   ├── nvt.mdp                 # NVT equilibration
│   ├── npt.mdp                 # NPT equilibration
│   ├── smd.mdp                 # Steered MD parameters
│   └── umbrella.mdp            # Umbrella sampling parameters
└── alignment/                  # Pulling direction results
    └── aligned_complex.pdb     # Aligned structure for visualization
```

---

## Step 2: Run Steered MD

After system building completes:

```bash
cd CDK2-CVT313/GMX_PROLIG_MD
bash run_smd.sh
```

This performs:

1. Energy minimization
2. NVT equilibration (heating to 300 K)
3. NPT equilibration (pressure coupling)
4. Steered MD: pulls the ligand out of the ATP-binding pocket along +Z at constant velocity

The SMD trajectory shows the ligand progressively leaving the binding site.

---

## Step 3: Run Umbrella Sampling

```bash
bash run_umbrella.sh
```

This extracts configurations from the SMD trajectory at regular intervals (default 0.12 nm spacing) and runs restrained simulations at each window.

For cluster submission, umbrella windows are independent and can be parallelized:

```bash
#!/bin/bash
#SBATCH --job-name=CDK2_umbrella
#SBATCH --array=0-39
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00

module load gromacs/2024.3

cd CDK2-CVT313/GMX_PROLIG_MD/umbrella/window_${SLURM_ARRAY_TASK_ID}
gmx mdrun -deffnm umbrella -v -ntmpi 1 -ntomp 8 -nb gpu -bonded gpu -pme gpu
```

---

## Step 4: WHAM Analysis

After all umbrella windows complete:

```bash
bash run_wham.sh
```

This runs the Weighted Histogram Analysis Method to reconstruct the PMF profile from the biased simulations. The output includes:

- **PMF profile**: Free energy as a function of the reaction coordinate (COM distance)
- **Histogram overlap**: Verifies sufficient sampling between adjacent windows
- **Bootstrap errors**: Statistical uncertainty from resampling

---

## Checking Results

### Verify Histogram Overlap

```bash
gmx wham -hist -f pullf-files.dat -it tpr-files.dat
xmgrace histo.xvg
```

Adjacent histograms should overlap significantly. Gaps indicate insufficient sampling - decrease `--umbrella-spacing` or increase `--umbrella-time`.

### Inspect the PMF Profile

```bash
xmgrace profile.xvg
```

A well-converged PMF should show:

- A **deep minimum** at the bound state (short COM distance)
- A **smooth rise** as the ligand leaves the pocket
- A **plateau** at large distances (bulk solvent)

The binding free energy is the difference between the plateau and the minimum:

$$
\Delta G_{\text{bind}} = W(\xi_{\text{bulk}}) - W(\xi_{\text{bound}})
$$

---

## Customization

### Adjust Umbrella Sampling Parameters

```bash
# Higher accuracy: closer windows, longer sampling
prism protein.pdb ligand.mol2 -ff amber14sb -lff gaff2 \
    --gaussian hf --isopt false -o CDK2-CVT313 --pmf \
    --umbrella-spacing 0.08 --umbrella-time 20

# Quick test: wider windows, shorter sampling
prism protein.pdb ligand.mol2 -ff amber14sb -lff gaff2 \
    --gaussian hf --isopt false -o CDK2-CVT313 --pmf \
    --umbrella-spacing 0.2 --umbrella-time 5
```

### Manual Pulling Direction

If you want to override the automatic direction optimization:

```bash
prism protein.pdb ligand.mol2 -ff amber14sb -lff gaff2 \
    --gaussian hf --isopt false -o CDK2-CVT313 --pmf \
    --pull-vector 100 200
```

Here `100` is the protein atom index near the pocket entrance and `200` is a ligand atom index.

### Use DFT Instead of HF

For potentially more accurate charges (at higher computational cost):

```bash
prism protein.pdb ligand.mol2 -ff amber14sb -lff gaff2 \
    --gaussian dft --isopt false -o CDK2-CVT313 --pmf
```

---

## Computational Cost

| Stage | Typical Time | Notes |
| --- | --- | --- |
| Gaussian RESP charges | 30 min - 2 hours | Depends on molecule size and basis set |
| System building + alignment | ~10 minutes | Includes MH optimization (~50k steps) |
| Steered MD | 4-24 hours | Single GPU |
| Umbrella sampling (40 windows) | 50-500 hours total | Embarrassingly parallel |
| WHAM analysis | < 1 hour | CPU only |

---

<div class="whats-next" markdown>

## What's Next

- [Read the PMF theory and options reference](../user-guide/pmf-calculations.md)
- [Submit umbrella windows to a cluster](../tutorials/cluster-submission.md)
- [Compare with MM/PBSA binding energy](../user-guide/analysis-tools.md#mmpbsa-binding-energy)

</div>
