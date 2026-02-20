# PMF Calculations

PRISM automates **Potential of Mean Force (PMF)** calculations for estimating protein-ligand binding free energies using steered MD, umbrella sampling, and WHAM analysis.

!!! example "Quick Start"
    ```bash
    prism protein.pdb ligand.mol2 -o output --pmf
    ```

## Overview

!!! info "What is PMF?"
    The Potential of Mean Force describes the free energy change along a reaction coordinate. For protein-ligand binding, PRISM pulls the ligand away from the binding pocket and calculates the binding free energy from the resulting energy profile.

The PRISM PMF workflow consists of three stages:

```mermaid
graph LR
    A[Equilibrated System] --> B[Steered MD]
    B --> C[Umbrella Sampling]
    C --> D[WHAM Analysis]
    D --> E[Binding Free Energy]
```

1. **Steered MD (SMD)**: Pull the ligand out of the binding pocket along an optimized direction
2. **Umbrella Sampling**: Restrain the ligand at evenly spaced windows along the unbinding pathway
3. **WHAM Analysis**: Reconstruct the PMF profile and extract binding free energy

## Step-by-Step Walkthrough

### Step 1: Build the PMF-Ready System

The `--pmf` flag tells PRISM to extend the simulation box along the pulling axis, align the system for optimal pulling, and generate SMD/umbrella MDP files:

```bash
prism protein.pdb ligand.mol2 -o pmf_output --pmf
```

PRISM automatically:

- Detects the binding pocket and selects an optimal pulling direction
- Aligns the system so the pulling axis is along the Z-axis
- Extends the box in the pulling direction to accommodate ligand unbinding
- Generates steered MD and umbrella sampling parameter files
- Creates run scripts for all stages

### Step 2: Run Steered MD

```bash
cd pmf_output/GMX_PROLIG_MD
bash run_smd.sh
```

This pulls the ligand from the bound state to a fully dissociated state.

### Step 3: Run Umbrella Sampling

After SMD completes, PRISM extracts configurations at regular intervals and sets up restrained simulations:

```bash
bash run_umbrella.sh
```

Each umbrella window runs independently and can be parallelized on a cluster.

### Step 4: Run WHAM Analysis

```bash
bash run_wham.sh
```

This produces the final PMF profile and binding free energy estimate.

## Options Reference

| Flag | Description | Default |
| --- | --- | --- |
| `--pmf` | Enable PMF calculation mode | off |
| `--pull-vector PROT LIG` | Atom indices defining the pulling direction (protein atom, ligand atom) | auto-detected |
| `--box-extension X Y Z` | Extra box length in each dimension (nm) | auto |
| `--umbrella-time` | Production time per umbrella window (ns) | `10.0` |
| `--umbrella-spacing` | Distance between umbrella windows (nm) | `0.12` |
| `--wham-begin` | Discard initial frames before WHAM (ps) | `1000` |
| `--wham-bootstrap` | Number of bootstrap samples for error estimation | `200` |

## Customizing the PMF Workflow

### Specifying the Pulling Direction

By default, PRISM uses a Metropolis-Hastings optimization to find the pulling direction with the fewest steric clashes. To override this with specific atom indices:

```bash
prism protein.pdb ligand.mol2 -o output --pmf --pull-vector 100 200
```

Here `100` is the index of a protein atom near the pocket entrance and `200` is a ligand atom index.

### Adjusting Umbrella Sampling

For higher accuracy, decrease the window spacing and increase the sampling time:

```bash
prism protein.pdb ligand.mol2 -o output --pmf \
  --umbrella-spacing 0.08 \
  --umbrella-time 20
```

For a quick test, use wider spacing and shorter sampling:

```bash
prism protein.pdb ligand.mol2 -o output --pmf \
  --umbrella-spacing 0.2 \
  --umbrella-time 5
```

### Extending the Box

If the ligand needs to be pulled further from the protein, extend the box explicitly:

```bash
prism protein.pdb ligand.mol2 -o output --pmf \
  --box-extension 0.0 0.0 5.0
```

## Cluster Submission

Umbrella windows are independent and can be submitted as an array job:

```bash
#!/bin/bash
#SBATCH --job-name=pmf_umbrella
#SBATCH --array=0-39
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00

module load gromacs/2024.3

cd pmf_output/GMX_PROLIG_MD/umbrella/window_${SLURM_ARRAY_TASK_ID}
gmx mdrun -deffnm umbrella -v -nb gpu -pme gpu
```

## Computational Resources

| Stage | Typical Time | Parallelization |
| --- | --- | --- |
| SMD | 4-24 hours (single GPU) | Single run |
| Umbrella | 50-500 hours total | Parallel windows |
| WHAM | < 1 hour (CPU) | Single run |

## Checking Convergence

After WHAM analysis, verify the results:

1. **Check histogram overlap:** Adjacent windows must overlap in configuration space
    ```bash
    gmx wham -hist -f pullf-files.dat -it tpr-files.dat
    ```

2. **Inspect the PMF profile:** It should be smooth without large jumps

3. **Check bootstrap errors:** The error from `--wham-bootstrap` should be < 1 kcal/mol

## Common Issues

!!! warning "Ligand pulls through protein"
    The pulling direction passes through the protein body instead of exiting the pocket.
    **Fix:** Use `--pull-vector` to specify the correct unbinding path, or let PRISM auto-detect it (the default uses Metropolis-Hastings optimization).

!!! warning "Poor WHAM convergence"
    Large gaps in the PMF profile or very high error bars.
    **Fix:** Decrease `--umbrella-spacing` and/or increase `--umbrella-time`.

!!! warning "Unrealistic binding energy"
    The computed binding energy does not match expectations.
    **Fix:** Ensure the system was properly equilibrated before PMF. Verify ligand parameterization. Try running multiple independent calculations.

## References

1. Torrie, G. M., & Valleau, J. P. (1977). *J. Comput. Phys.*, 23(2), 187-199.
2. Kumar, S., et al. (1992). *J. Comput. Chem.*, 13(8), 1011-1021.
3. Hub, J. S., et al. (2010). *J. Chem. Theory Comput.*, 6(10), 3713-3720.

<div class="whats-next" markdown>

## What's Next

- [Follow the PMF Tutorial for a step-by-step example](../tutorials/pmf-tutorial.md)
- [Learn about Running Simulations for equilibration before PMF](running-simulations.md)
- [See Analysis Tools for examining PMF trajectories](analysis-tools.md)

</div>
