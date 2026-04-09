# FEP Calculations

PRISM automates **Free Energy Perturbation (FEP)** calculations for relative binding free energies between similar ligands using GROMACS hybrid topologies.

!!! example "Quick Start"
    ```bash
    prism protein.pdb ligand_ref.mol2 -o fep_output \
      --fep \
      --mutant ligand_mut.mol2 \
      --ligand-forcefield gaff2 \
      --forcefield amber14sb_OL15 \
      --fep-config fep.yaml
    ```

If you want the full YAML parameter reference, jump directly to [FEP YAML Reference](#fep-yaml-reference).

## Overview

PRISM's FEP workflow consists of five stages:

```mermaid
graph LR
    A[Ligand Pair] --> B[Atom Mapping]
    B --> C[Hybrid Topology]
    C --> D[FEP Scaffold]
    D --> E[Bound Leg]
    D --> F[Unbound Leg]
    E --> G[Analysis]
    F --> G
```

1. **Atom Mapping**: distance-based matching between the reference and mutant ligands
2. **Hybrid Topology**: one ligand topology with A/B-state parameters
3. **FEP Scaffold**: bound and unbound system setup, MDP generation, and run scripts
4. **Lambda Windows**: per-window equilibration and production inputs
5. **Analysis**: BAR, MBAR, and TI estimators over the window data

## Mapping and Hybrid Topology

### Distance-Based Atom Mapping

PRISM maps atoms using geometry and element identity, with optional charge-based filtering.

Default mapping parameters:

| Parameter | Default | Meaning |
|---|---:|---|
| `dist_cutoff` | `0.6 nm` | Distance cutoff for candidate atom matches. |
| `charge_cutoff` | `0.05 e` | Charge-difference filter used during mapping/classification. |
| `charge_common` | `mean` | Default shared-atom charge strategy. |
| `charge_reception` | `surround` | Default redistribution target for compensating charge. |

For detailed parameter descriptions, see the [`mapping` section](#mapping-section) in the FEP YAML Reference.

!!! note
    PRISM uses GROMACS-style length units internally. A mapping cutoff of `0.6` means **0.6 nm**, not 0.6 Å.

### Atom Classes

PRISM classifies atoms as:

- **Common atoms**: shared between both ligands
- **Transformed atoms**: only present in state A or only in state B
- **Surrounding atoms**: retained in the hybrid representation but handled with state-specific parameters near the mutation region

### Charge Redistribution

**Why charge redistribution matters**: FEP calculations measure free energy differences by alchemically transforming one ligand into another. In practice, smaller and more local perturbations usually improve phase-space overlap and make the calculation easier to converge.

When atoms are classified as **common** (shared between both ligands), keeping their electrostatic properties as consistent as possible usually reduces the magnitude of the alchemical transformation. This often leads to:

- **Better convergence** — fewer lambda windows needed for adequate overlap
- **Lower variance** — smaller statistical uncertainty in the final ΔG
- **Improved physical realism** — common regions should behave consistently in both states

PRISM implements this through two complementary mechanisms:

1. **`charge_common`**: Controls how shared atoms inherit charge information
2. **`charge_reception`**: When atoms transform, PRISM redistributes the compensating charge so that the end states remain consistent with the intended charge model

The following figure visualizes these strategies:

![Charge redistribution modes](../assets/fep/charge_redistribution.png)

For parameter details and available options, see the [`mapping` section](#mapping-section) below.

For the current GROMACS discussion of free-energy pathways, end states, and interaction handling, see the current GROMACS free-energy implementation and interaction references.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

### Hybrid Topology

The generated `hybrid.itp` stores both state A and state B parameters in a single molecule description. In GROMACS this means the **end states** are encoded in the topology, while the **path** between them is controlled by the lambda schedule in the generated MDP files.

PRISM-generated MDPs write separate lambda vectors for:

| Lambda vector | Role |
|---|---|
| `coul-lambdas` | Controls electrostatic transformation. |
| `vdw-lambdas` | Controls van der Waals transformation. |
| `bonded-lambdas` | Controls bonded-term interpolation. |
| `mass-lambdas` | Controls mass interpolation. |

So the implementation is more specific than a simple scalar interpolation formula.

For the current GROMACS description of free-energy end states, lambda vectors, and interaction handling, see the current manual and the free-energy reference sections.[^gmx-current] [^gmx-free-energy-impl] [^gmx-free-energy-interactions] The `gmx grompp` manual also documents `-rb` for B-state reference coordinates.[^gmx-grompp]

## Lambda Schedules

PRISM currently supports three lambda schedule strategies. The exact keys and defaults are listed in the [`lambda` section](#lambda-section) of the FEP YAML Reference.

### Supported schedule modes

| Mode | Default? | Meaning |
|---|---|---|
| `decoupled` | yes | Coulomb transforms first, then VDW transforms. |
| `coupled` | no | Coulomb and VDW follow the same schedule together from 0 to 1. |
| `custom` | no | User supplies explicit `custom_coul_lambdas` and `custom_vdw_lambdas`. Shorter arrays are padded with their final value. |

The following figure visualizes the three lambda schedule modes:

![PRISM lambda schedule modes](../assets/fep/lambda_schedule_modes.png)

**Figure explanation**:
- **Decoupled (default)**: Shows the two-stage transformation where Coulomb (red) completes first (windows 0-11), followed by VDW (blue) (windows 12-31). Bonded (gray dashed) and mass (green dotted) parameters follow the Coulomb schedule.
- **Coupled**: Coulomb and VDW transform simultaneously across all 32 windows, maintaining identical lambda values throughout.
- **Custom**: Example showing user-defined arrays with different transition points — Coulomb jumps early (window 2), VDW transitions later (window 3).

### Supported distributions

| Distribution | Default? | Meaning |
|---|---|---|
| `linear` | no | Evenly spaced lambda points. |
| `nonlinear` | yes | Denser sampling near the endpoints. |
| `quadratic` | no | Endpoint-biased spacing based on a quadratic shape. |

Current GROMACS references for lambda schedules and interaction-specific lambda control are the free-energy implementation and interaction chapters.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

## Build Workflow

### Build the Scaffold

```bash
prism protein.pdb ligand_ref.mol2 -o fep_output \
  --fep \
  --mutant ligand_mut.mol2 \
  --ligand-forcefield gaff2 \
  --forcefield amber14sb_OL15 \
  --fep-config fep.yaml
```

PRISM then:

1. maps the ligand pair
2. builds the hybrid topology and `hybrid.gro`
3. prepares the **bound** leg (protein + hybrid ligand)
4. prepares the **unbound** leg (hybrid ligand in solvent)
5. writes per-window MDPs
6. writes a root-level `run_fep.sh`
7. writes per-leg helpers such as `run_prod_standard.sh` and `run_prod_repex.sh`
8. writes `fep_scaffold.json`

### Generated Window Workflow

For each leg, PRISM currently prepares:

1. leg-level EM
2. leg-level NVT
3. leg-level NPT
4. per-window `em_short`
5. per-window `npt_short`
6. per-window production (`prod`)

This means each lambda window is not launched directly from the leg-level NPT state; it first receives its own short lambda-specific relaxation.

This is PRISM's generated workflow for robustness and restartability. It is an implementation choice, not a GROMACS-imposed universal FEP sequence.

### Generated MDP Behavior

The generated per-window MDPs currently include:

| MDP item | Current behavior |
|---|---|
| `init-lambda-state` | Set to the current window index. |
| `calc-lambda-neighbors` | `-1` in window-level equilibration and production templates. |
| Lambda vectors | Explicit `coul-lambdas`, `vdw-lambdas`, `bonded-lambdas`, and `mass-lambdas` are written. |
| Soft-core settings | `sc-alpha` and `sc-sigma` are included in the generated templates. |

For the current GROMACS description of these free-energy controls, see the free-energy implementation and interaction chapters.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

## Running FEP

```bash
cd fep_output/GMX_PROLIG_FEP

bash run_fep.sh bound
bash run_fep.sh unbound
bash run_fep.sh all
```

If the scaffold contains multiple replicas, the root script also supports:

```bash
bash run_fep.sh bound1
bash run_fep.sh unbound2
bash run_fep.sh bound1-3
```

### Production Modes

| Mode | Meaning | Notes |
|---|---|---|
| `standard` | Run windows independently. | `parallel_windows` controls concurrency; GPUs are assigned round-robin across windows. |
| `repex` | Run lambda replica exchange. | PRISM writes `run_prod_repex.sh` and launches production with `gmx_mpi -multidir`. |

### Smoke Testing

For infrastructure checks only:

```bash
export PRISM_MDRUN_NSTEPS=100
bash run_fep.sh bound
bash run_fep.sh unbound
unset PRISM_MDRUN_NSTEPS
```

!!! warning
    `PRISM_MDRUN_NSTEPS` is for infrastructure validation only. Do not use these shortened runs for scientific free-energy estimates.

## Analysis

PRISM analyzes FEP results using three established free energy estimators:

| Estimator | Full Name | Description | Use Case |
|---|---|---|---|
| **TI** | Thermodynamic Integration | Numerical integration of ∂H/∂λ across windows | Good for smooth transformations with well-sampled gradients |
| **BAR** | Bennett Acceptance Ratio | Optimal overlap between neighboring windows | Robust for standard FEP with adequate sampling |
| **MBAR** | Multistate BAR | Uses all data simultaneously with reweighting | Most efficient; provides overlap matrix for quality control |

**How it works**:
1. PRISM reads `dhdl.xvg` files from each `window_*` directory (contains ∂H/∂λ time series)
2. For each leg (bound/unbound), computes ΔG and uncertainty using the selected estimator(s)
3. Binding free energy: ΔG_bind = ΔG_unbound - ΔG_bound
4. Bootstrap analysis provides error estimates (controlled by `--bootstrap-n-jobs`)
5. HTML report includes overlap matrices, convergence plots, and estimator comparison

**Backend options**:
- `alchemlyb` (default): Python library supporting TI/BAR/MBAR; required for multi-estimator mode
- `gmx_bar`: Native GROMACS `gmx bar` command; BAR-only, useful for validation

Example:

```bash
prism --fep-analyze \
  --bound-dir fep_output/GMX_PROLIG_FEP/bound/repeat1 \
  --unbound-dir fep_output/GMX_PROLIG_FEP/unbound/repeat1 \
  --estimator MBAR BAR TI \
  --bootstrap-n-jobs 8 \
  --output fep_results.html \
  --json fep_results.json
```

## Practical Notes

- Prefer an explicit `fep.yaml` for reproducible lambda schedules.
- The CLI still exposes `--lambda-windows`, but it is a simplified override. If you care about exact schedules, define the lambda section explicitly in YAML.
- `execution.mode`, `parallel_windows`, `num_gpus`, `total_cpus`, and `omp_threads` only affect generated run scripts; they do not change the topology or lambda schedule.
- If you use `custom` lambda schedules, validate the arrays carefully before large production campaigns.

## FEP YAML Reference

If you started from the workflow sections above and now want the exact configuration keys, use the following reference.

### Configuration Map

```mermaid
block-beta
    columns 3

    B["mapping<br/>dist_cutoff<br/>charge_cutoff<br/>charge_common<br/>charge_reception"]
    C["lambda<br/>strategy / distribution<br/>windows / coul_windows / vdw_windows<br/>custom_coul_lambdas / custom_vdw_lambdas"]
    D["simulation<br/>equilibration / production time<br/>dt / temperature / pressure<br/>soft_core: alpha / sigma"]

    F["electrostatics<br/>rcoulomb"]
    A(("fep.yaml"))
    G["vdw<br/>rvdw"]

    H["output<br/>trajectory / energy / log intervals"]
    I["execution<br/>mode / num_gpus / parallel_windows<br/>total_cpus / omp_threads<br/>use_gpu_pme / mdrun_update_mode"]
    J["replicas"]

    A --> B
    A --> C
    A --> D
    A --> F
    A --> G
    A --> H
    A --> I
    A --> J

    style A fill:#e8f1ff,stroke:#4c78a8,stroke-width:2px,color:#111
    style B fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style C fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style D fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style F fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style G fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style H fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style I fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
    style J fill:#f8f9fb,stroke:#7a869a,stroke-width:1px,color:#111
```

Mermaid flowcharts do not provide a true radial layout. This map therefore uses Mermaid's block diagram syntax to keep a center node with sections above and below while preserving outward arrows.[^mermaid-block]

### Example `fep.yaml`

```yaml
mapping:
  dist_cutoff: 0.6
  charge_cutoff: 0.05
  charge_common: mean
  charge_reception: surround

lambda:
  strategy: decoupled
  distribution: nonlinear
  windows: 32
  coul_windows: 12
  vdw_windows: 20
  # custom_coul_lambdas: [0.0, 0.2, 0.5, 1.0]
  # custom_vdw_lambdas: [0.0, 0.0, 0.5, 1.0]

simulation:
  equilibration_nvt_time_ps: 500
  equilibration_npt_time_ps: 500
  production_time_ns: 5.0
  dt: 0.002
  temperature: 310
  pressure: 1.0

soft_core:
  alpha: 0.5
  sigma: 0.3

electrostatics:
  rcoulomb: 1.0

vdw:
  rvdw: 1.0

output:
  trajectory_interval_ps: 500
  energy_interval_ps: 10
  log_interval_ps: 10

execution:
  mode: standard
  num_gpus: 4
  parallel_windows: 4
  total_cpus: 56
  omp_threads: 14
  use_gpu_pme: true
  mdrun_update_mode: auto

replicas: 1
```

### `mapping` section

| Key | Default | Meaning |
|---|---:|---|
| `dist_cutoff` | `0.6` | Distance cutoff for atom mapping, in **nm**. |
| `charge_cutoff` | `0.05` | Charge-difference cutoff used during mapping/classification. |
| `charge_common` | `mean` | How shared atoms inherit charge information: `ref`, `mut`, `mean`. |
| `charge_reception` | `surround` | Where compensating charge is redistributed after common-atom assignment. Common default is `surround`. |

#### Charge redistribution modes

Current user-facing workflows currently validate these `charge_reception` modes:

| Mode | Default? | Meaning |
|---|---|---|
| `surround` | yes | Redistribute onto atoms surrounding the transformed region. |
| `unique` | no | Redistribute only onto atoms unique to the transformed region. |
| `none` | no | Disable redistribution after common-atom assignment. |

!!! note
    This figure is also used in the PRISM preprint and is included here as a conceptual aid. Exact atom membership in each redistribution set still depends on the current topology and mapping result.

!!! warning
    The lower-level hybrid-topology code contains an experimental `surround_ext` redistribution branch, but the standard distance-based mapping path currently validates only `surround`, `unique`, and `none`. Treat `surround_ext` as internal/experimental rather than a normal documented workflow option.

### `lambda` section

| Key | Default | Meaning |
|---|---:|---|
| `strategy` | `decoupled` | Lambda scheduling mode: `decoupled`, `coupled`, or `custom`. |
| `distribution` | `nonlinear` | Distribution along the schedule: `linear`, `nonlinear`, or `quadratic`. |
| `windows` | `32` | Total window count for `coupled`, and target total for `decoupled`. |
| `coul_windows` | `12` | Coulomb-stage window count in `decoupled` mode. |
| `vdw_windows` | `20` | VDW-stage window count in `decoupled` mode. |
| `custom_coul_lambdas` | none | Required for `custom`; explicit Coulomb lambda array. |
| `custom_vdw_lambdas` | none | Required for `custom`; explicit VDW lambda array. |

If you use only the top-level CLI flag `--lambda-windows`, treat it as a simple convenience override. For reproducible FEP work, prefer an explicit YAML lambda section.

Reference: see the current GROMACS free-energy implementation and interaction chapters.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

### `simulation` section

| Key | Default | Meaning |
|---|---:|---|
| `equilibration_nvt_time_ps` | `500` | Leg-level NVT equilibration time in ps. |
| `equilibration_npt_time_ps` | `500` | Leg-level NPT equilibration time in ps. |
| `production_time_ns` | `5.0` | Production length per lambda window in ns. |
| `dt` | `0.002` | MD time step in ps. |
| `temperature` | `310` | Simulation temperature in K. |
| `pressure` | `1.0` | Simulation pressure in bar. |

### `soft_core` section

| Key | Default | Meaning |
|---|---:|---|
| `alpha` | `0.5` | Soft-core alpha parameter for alchemical nonbonded interactions. |
| `sigma` | `0.3` | Soft-core sigma parameter in nm. |

### `electrostatics` and `vdw` sections

| Section | Key | Default | Meaning |
|---|---|---:|---|
| `electrostatics` | `rcoulomb` | `1.0` | Requested Coulomb cutoff in nm. |
| `vdw` | `rvdw` | `1.0` | Requested VDW cutoff in nm. |

!!! note
    For hybrid FEP MDP generation, the runtime templates clamp `rcoulomb` and `rvdw` to at least `1.4 nm` for stability even if your YAML requests smaller values.

### `output` section

| Key | Default | Meaning |
|---|---:|---|
| `trajectory_interval_ps` | `500` | Trajectory output interval in ps. |
| `energy_interval_ps` | `10` | Energy output interval in ps. |
| `log_interval_ps` | `10` | Log output interval in ps. |

### `execution` section

| Key | Default | Meaning |
|---|---:|---|
| `mode` | `standard` | Production execution mode: `standard` or `repex`. |
| `num_gpus` | `null` | Total GPUs available to generated scripts. |
| `parallel_windows` | `null` | Concurrent windows in `standard` mode; if omitted, scripts fall back to GPU count. |
| `total_cpus` | `null` | Total CPUs available; used to derive OpenMP threads when `omp_threads` is not set. |
| `omp_threads` | `null` | Manual OpenMP thread override per worker/rank. |
| `use_gpu_pme` | `true` | Whether production helpers request `-pme gpu`. |
| `mdrun_update_mode` | `auto` | Runtime update mode: `auto`, `gpu`, `cpu`, or `none`. |
| `use_gpu_update` | `false` | Legacy compatibility switch; current FEP defaults keep GPU-side update off unless explicitly enabled. |

### Top-level key

| Key | Default | Meaning |
|---|---:|---|
| `replicas` | `1` | Number of bound/unbound repeat directories to scaffold. |

## References

[^gmx-current]: GROMACS current manual. https://manual.gromacs.org/current/
[^gmx-free-energy-impl]: GROMACS current reference manual, *Free energy implementation*. https://manual.gromacs.org/current/reference-manual/special/free-energy-implementation.html
[^gmx-free-energy-interactions]: GROMACS current reference manual, *Free-energy interactions*. https://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html
[^gmx-grompp]: GROMACS current online help, `gmx grompp`. https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html
[^mermaid-block]: Mermaid block diagrams. https://mermaid.js.org/syntax/block.html

## Troubleshooting

### Mapping looks chemically wrong

- verify protonation and tautomer states first
- prefer MOL2 or SDF when PDB loses bond-order information
- inspect `common/hybrid/mapping.html` before trusting a scaffold

### Windows fail early

- inspect leg-level `build/em.log`, `build/nvt.log`, and `build/npt.log`
- inspect the hybrid files in `common/hybrid/`
- rerun a smoke test before launching long production

### Analysis finds no windows

- pass the repeat directory, not only the FEP root
- confirm that each `window_*` directory contains production outputs such as `prod.*` and `dhdl.xvg`

## See Also

- [FEP Tutorial](../tutorials/fep-tutorial.md)
- [Analysis Tools](analysis-tools.md)
- [Configuration](configuration.md)
