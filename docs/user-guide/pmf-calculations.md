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

## Theory

### Pulling Direction Optimization

A critical challenge in PMF calculations is selecting a pulling direction that avoids steric clashes with the protein. PRISM solves this by formulating the direction selection as an optimization problem on the unit sphere $S^2$.

For a binding pocket with $N$ protein atoms at positions $\{p_i\}$ and a ligand with $M$ atoms at positions $\{l_j\}$, the optimal pulling direction $\mathbf{v}^*$ maximizes the clearance between the pulling path and all pocket atoms:

$$
\mathbf{v}^* = \arg\min_{\mathbf{v} \in S^2} \left[ -\sum_{i=1}^{N} \min_{j=1}^{M} \| (p_i - l_j) \times \mathbf{v} \| \right]
$$

Here $\| (p_i - l_j) \times \mathbf{v} \|$ is the perpendicular distance from pocket atom $p_i$ to the line passing through ligand atom $l_j$ in direction $\mathbf{v}$. By summing the minimum such distance over all pocket atoms, the objective function penalizes directions where pulling lines pass close to any pocket atom.

This non-convex optimization is solved using **Metropolis-Hastings simulated annealing** on the unit sphere:

1. Initialize with the protein COM → ligand COM direction
2. At each step, propose a random perturbation on $S^2$ (spherical coordinates with step size $\sigma = 0.15$ rad)
3. Accept or reject based on the Metropolis criterion with a cooling schedule: $T_k = T_0 \cdot r^k$ where $T_0 = 2.0$ and $r = 0.99995$
4. Converge when the relative standard deviation of the energy over the last 8000 iterations falls below $5 \times 10^{-4}$

After optimization, the complex is rotated so that $\mathbf{v}^*$ aligns with the +Z axis, enabling standard GROMACS pulling along Z.

### Steered Molecular Dynamics (SMD)

In steered MD, an external harmonic potential is applied to the ligand center of mass and moved at constant velocity along the pulling direction:

$$
U_{\text{pull}}(t) = \frac{1}{2} k \left[ z_{\text{lig}}(t) - z_0 - v \cdot t \right]^2
$$

where $k$ is the spring constant (typically 1000 kJ/mol/nm²), $z_{\text{lig}}(t)$ is the ligand COM position, $z_0$ is the initial position, and $v$ is the pull rate. The work performed on the system is recorded and used to generate starting configurations for umbrella sampling.

### Umbrella Sampling and WHAM

Direct computation of the PMF from a single SMD trajectory is unreliable due to non-equilibrium effects. Instead, **umbrella sampling** places harmonic restraints at evenly spaced windows $\{\xi_i\}$ along the reaction coordinate:

$$
U_i(\xi) = \frac{1}{2} k_i (\xi - \xi_i)^2
$$

Each window is simulated independently, producing a biased probability distribution $P_i^{\text{biased}}(\xi)$.

The **Weighted Histogram Analysis Method (WHAM)** reconstructs the unbiased PMF $W(\xi)$ by solving the coupled self-consistency equations:

$$
W(\xi) = -k_B T \ln \left[ \frac{\sum_{i=1}^{S} n_i \, h_i(\xi)}{\sum_{i=1}^{S} n_i \, \exp\!\left[-\beta\!\left(U_i(\xi) - F_i\right)\right]} \right]
$$

$$
\exp(-\beta F_i) = \int \exp\!\left[-\beta\!\left(W(\xi) + U_i(\xi)\right)\right] d\xi
$$

where $S$ is the number of windows, $n_i$ is the number of samples in window $i$, $h_i(\xi)$ is the histogram count, $\beta = 1/k_B T$, and $F_i$ is the free energy constant for each window. These equations are solved iteratively until convergence, and the binding free energy is:

$$
\Delta G_{\text{bind}} = W(\xi_{\text{bound}}) - W(\xi_{\text{bulk}})
$$

Statistical uncertainties are estimated via Bayesian bootstrap resampling of the histograms.

---

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
