# FEP Workflow Tutorial

**Level:** Advanced  
**Time:** 2-4 hours (+ computation time)  
**Topics:** Relative binding free energy, hybrid topologies, lambda schedules, execution modes, analysis

This tutorial walks through the current PRISM FEP workflow using the **42-38** HIF-2α ligand pair.

## Objectives

You will learn how to:

- prepare a ligand pair for PRISM FEP
- inspect atom mapping and the hybrid topology
- understand how PRISM transforms the system across lambda windows
- run the generated bound and unbound legs
- choose between `standard` and `repex` execution
- analyze completed window data with BAR, MBAR, and TI

## Prerequisites

- completed [Basic Tutorial](basic-tutorial.md) or equivalent PRISM familiarity
- GROMACS installed and available as `gmx`
- enough CPU/GPU resources for your chosen window schedule
- basic understanding of alchemical free-energy calculations

!!! note
    Multiple GPUs improve throughput, but PRISM can still generate valid scaffolds for smaller GPU counts or CPU-only debugging workflows.

## Background: The 42-38 System

We use the **42-38** ligand pair from HIF-2α:

- **Ligand 42**: reference ligand (state A)
- **Ligand 38**: mutant ligand (state B)

This pair is suitable for FEP because the chemical change is local and the ligands share a clear common scaffold.

## Step 1: Obtain the Example System

```bash
git clone https://github.com/AIB001/PRISM.git
cd PRISM/tests/gxf/FEP/unit_test/42-38
```

Relevant files:

```text
42-38/
├── input/
│   ├── sys.pdb
│   ├── 42.pdb
│   └── 38.pdb
└── configs/
    ├── fep_gaff2.yaml
    └── case_gaff2.yaml
```

## Step 2: Review the FEP Configuration

```bash
cat configs/fep_gaff2.yaml
```

A representative setup is:

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

simulation:
  equilibration_nvt_time_ps: 500
  equilibration_npt_time_ps: 500
  production_time_ns: 2.0
  dt: 0.002
  temperature: 310
  pressure: 1.0

execution:
  mode: standard
  num_gpus: 4
  parallel_windows: 4
  total_cpus: 56
  use_gpu_pme: true
  mdrun_update_mode: auto

replicas: 1
```

Key points:

| Item | Current behavior |
|---|---|
| `dist_cutoff: 0.6` | Means **0.6 nm**, not 0.6 Å. |
| `strategy: decoupled` | This is the current code default. |
| Default windows | `32` total windows split as `12` Coulomb + `20` VDW windows. |
| CLI vs YAML | Use an explicit FEP YAML when you care about exact schedules. |

## Step 3: Understand What the Lambda Windows Actually Do

This is the most important conceptual point.

In PRISM's default **decoupled** schedule:

1. **Coulomb stage**
   - `coul-lambdas` move from `0 -> 1`
   - `vdw-lambdas` stay at `0`
2. **VDW stage**
   - `coul-lambdas` stay at `1`
   - `vdw-lambdas` move from `0 -> 1`

At the same time:

| Lambda vector | Current behavior |
|---|---|
| `bonded-lambdas` | Follow the Coulomb schedule. |
| `mass-lambdas` | Follow the Coulomb schedule. |

So PRISM does **not** simply scale a single scalar lambda through the whole transformation. It writes separate lambda vectors into each MDP.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

If you choose:

- `strategy: coupled`
  - Coulomb and VDW transform together across all windows
- `strategy: custom`
  - you must provide both `custom_coul_lambdas` and `custom_vdw_lambdas`
  - if one custom array is shorter than the other, PRISM pads the shorter one with its final value

The following figure visualizes the three lambda schedule strategies:

![PRISM lambda schedule modes](../assets/fep/lambda_schedule_modes.png)

**What the figure shows**:
- **Decoupled**: The default two-stage approach — electrostatics (red line) finish first (windows 0-11), then van der Waals (blue squares) take over (windows 12-31)
- **Coupled**: Both red and blue lines move together from 0 to 1 across all windows
- **Custom**: A user-defined example where you control exactly when each term transforms

Supported schedule distributions are:

| Distribution | Meaning |
|---|---|
| `linear` | Evenly spaced lambda points. |
| `nonlinear` | Denser endpoint sampling; current default. |
| `quadratic` | Endpoint-biased spacing from a quadratic curve. |

### Interpreting charge redistribution

**The core principle**: In practice, smaller and more local perturbations usually improve phase-space overlap and make FEP easier to converge. By keeping electrostatic properties consistent in the common (shared) regions of the ligands, we often:

- **Reduce the magnitude of the perturbation** → better convergence
- **Lower statistical uncertainty** → more reliable ΔG estimates
- **Maintain physical consistency** → common regions behave similarly in both states

The following figure visualizes the charge redistribution strategies:

![Charge redistribution modes](../assets/fep/charge_redistribution.png)

For parameter details, see the table in [Which options should users know about first?](#which-options-should-users-know-about-first).

## Step 4: Build the FEP Scaffold

```bash
prism input/sys.pdb input/42.pdb -o amber14sb_OL15-mut_gaff2 \
  --fep \
  --mutant input/38.pdb \
  --ligand-forcefield gaff2 \
  --forcefield amber14sb_OL15 \
  --fep-config configs/fep_gaff2.yaml \
  --config configs/case_gaff2.yaml
```

Expected outputs:

- `amber14sb_OL15-mut_gaff2/GMX_PROLIG_FEP/`
- `common/hybrid/` with hybrid topology and mapping artifacts
- `bound/repeat1/` and `unbound/repeat1/`
- root-level `run_fep.sh`
- `fep_scaffold.json`

## Step 5: Inspect Atom Mapping and Hybrid Files

Open the mapping report:

```bash
firefox amber14sb_OL15-mut_gaff2/GMX_PROLIG_FEP/common/hybrid/mapping.html
```

Inspect the scaffold summary:

```bash
python -m json.tool amber14sb_OL15-mut_gaff2/GMX_PROLIG_FEP/fep_scaffold.json | less
```

Check that:

1. the transformed region is local and chemically sensible
2. the common scaffold matches your expectation
3. there are no obviously incorrect long-range swaps
4. the manifest contains a non-empty `mapping` section

## Step 6: Inspect the Scaffold Layout

```bash
cd amber14sb_OL15-mut_gaff2/GMX_PROLIG_FEP
find . -maxdepth 3 | sed -n '1,80p'
```

Typical layout:

```text
GMX_PROLIG_FEP/
├── run_fep.sh
├── fep_scaffold.json
├── common/
│   └── hybrid/
├── bound/
│   └── repeat1/
│       ├── build/
│       ├── input/
│       ├── mdps/
│       ├── run_prod_standard.sh
│       ├── run_prod_repex.sh
│       └── window_00/ ... window_31/
└── unbound/
    └── repeat1/
        ├── build/
        ├── input/
        ├── mdps/
        ├── run_prod_standard.sh
        ├── run_prod_repex.sh
        └── window_00/ ... window_31/
```

Important points:

- `run_fep.sh` drives leg-level EM, NVT, NPT, then production windows
- `mdps/` contains leg-level and window-level MDP files
- `window_*` directories are the per-window production workspaces
- each window gets its own short lambda-specific relaxation before production

## Step 7: Run a Short Smoke Test

Before full production, validate the scaffold with a short test:

```bash
export PRISM_MDRUN_NSTEPS=100
bash run_fep.sh bound
bash run_fep.sh unbound
unset PRISM_MDRUN_NSTEPS
```

This is useful to validate:

| Check | Why it matters |
|---|---|
| Hybrid topology | Confirms the scaffold is chemically coherent. |
| Generated scripts | Confirms the execution layer was written correctly. |
| Per-window grompp/mdrun setup | Confirms lambda-specific files and commands are valid. |
| GPU and CPU allocation logic | Confirms runtime resource assignment behaves as expected. |

!!! warning
    `PRISM_MDRUN_NSTEPS` is only for infrastructure validation. Do not use these truncated runs for scientific free-energy estimation.

## Step 8: Choose an Execution Mode

PRISM currently supports two production modes.

### Standard mode

This is the usual default:

```yaml
execution:
  mode: standard
  num_gpus: 4
  parallel_windows: 4
```

Behavior:

| Property | `standard` mode behavior |
|---|---|
| Window execution | Windows are independent. |
| Concurrency | Up to `parallel_windows` windows run concurrently. |
| GPU assignment | GPUs are assigned round-robin. |
| Best for | Simple throughput-oriented execution. |

### Replica-exchange mode

```yaml
execution:
  mode: repex
  num_gpus: 4
```

Behavior:

| Property | `repex` mode behavior |
|---|---|
| Helper script | PRISM writes `run_prod_repex.sh`. |
| Launch style | Production is launched with `gmx_mpi -multidir`. |
| Best for | Lambda replica exchange instead of independent windows. |

!!! note
    `execution.parallel_windows` matters in `standard` mode. In `repex`, the main switch is `execution.mode: repex`.

## Step 9: Run Full Production

After the smoke test passes:

```bash
bash run_fep.sh bound
bash run_fep.sh unbound
```

Or both legs at once:

```bash
bash run_fep.sh all
```

If your scaffold contains multiple replicas, the master script also supports:

```bash
bash run_fep.sh bound1
bash run_fep.sh unbound2
bash run_fep.sh bound1-3
```

## Step 10: Monitor Progress

Examples:

```bash
find bound/repeat1 -maxdepth 1 -type d -name 'window_*' | wc -l
ls bound/repeat1/window_00/
```

What to check:

| Item | Expected sign |
|---|---|
| EM | Completes with a reasonable final force. |
| NVT/NPT | Produce `*.gro` and `*.cpt`. |
| Window outputs | Each `window_*` directory produces `prod.*` and `dhdl.xvg`. |
| Logs | Contain `Finished mdrun on rank 0`. |

## Step 11: Analyze Completed Windows

```bash
prism --fep-analyze \
  --bound-dir amber14sb_OL15-mut_gaff2/GMX_PROLIG_FEP/bound/repeat1 \
  --unbound-dir amber14sb_OL15-mut_gaff2/GMX_PROLIG_FEP/unbound/repeat1 \
  --estimator MBAR BAR TI \
  --bootstrap-n-jobs 8 \
  --output fep_results.html \
  --json fep_results.json
```

This reads `dhdl.xvg` files from each `window_*` directory and computes binding free energies using three estimators:

| Estimator | What it does | When to use |
|---|---|---|
| **TI** (Thermodynamic Integration) | Integrates ∂H/∂λ across windows | Smooth transformations with well-sampled gradients |
| **BAR** (Bennett Acceptance Ratio) | Compares neighboring window pairs | Standard FEP with adequate sampling |
| **MBAR** (Multistate BAR) | Reweights all data simultaneously | Most efficient; provides overlap matrix for quality checks |

**Analysis workflow**:
1. Reads `dhdl.xvg` time series from each window
2. Computes ΔG_bound and ΔG_unbound (with bootstrap error estimates)
3. Calculates binding free energy: ΔG_bind = ΔG_unbound - ΔG_bound
4. Generates HTML report with overlap matrices and convergence plots

## Step 12: Validate the Results

Open the report:

```bash
firefox fep_results.html
```

Check:

| Metric | What to look for |
|---|---|
| Overlap matrix | Neighboring windows overlap adequately. |
| Estimator agreement | BAR, MBAR, and TI are broadly consistent. |
| Convergence | Estimates stabilize as more data are accumulated. |
| Uncertainty | Bootstrap error bars are reasonable for the sampling length. |

## Common Questions

### Why are there so many window-level MDP files?

Because PRISM writes explicit per-window lambda schedules, plus short per-window equilibration stages (`em_short`, `npt_short`) before production.

### Is `--lambda-windows` enough by itself?

Not if you want reproducibility. For serious FEP work, define `lambda.strategy`, `distribution`, and any stage counts explicitly in `fep.yaml`.

### Why does PRISM redistribute charges during mapping?

FEP accuracy depends on minimizing the perturbation between ligands. If common (shared) atoms have very different charges in the two states, the alchemical transformation becomes larger and harder to converge:

- **Larger perturbation** → poorer convergence → higher uncertainty in ΔG
- **Smaller perturbation** → better convergence → more reliable results

By using `charge_common: mean` and `charge_reception: surround`, PRISM keeps electrostatic properties consistent in shared regions, confining the major changes to the actual mutation site.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

### Which options should users know about first?

These are the most important knobs:

| Parameter | Why it matters |
|---|---|
| `mapping.dist_cutoff` | Controls how permissive atom mapping is. |
| `mapping.charge_common` | Changes how shared atoms inherit charge information. |
| `mapping.charge_reception` | Changes where balancing charge is redistributed. |
| `lambda.strategy` | Chooses decoupled, coupled, or custom lambda behavior. |
| `lambda.distribution` | Controls where windows are denser or sparser. |
| `lambda.windows` | Sets the total schedule size. |
| `lambda.coul_windows` | Controls the Coulomb stage length in decoupled mode. |
| `lambda.vdw_windows` | Controls the VDW stage length in decoupled mode. |
| `simulation.production_time_ns` | Sets per-window production length. |
| `execution.mode` | Chooses standard versus replica-exchange execution. |
| `execution.parallel_windows` | Sets standard-mode concurrency. |
| `execution.num_gpus` | Affects GPU allocation in generated scripts. |
| `execution.total_cpus` / `execution.omp_threads` | Control CPU/OpenMP layout. |
| `execution.use_gpu_pme` | Requests GPU PME in production helpers. |
| `execution.mdrun_update_mode` | Controls GPU/CPU update behavior in runtime commands. |
| `replicas` | Controls how many repeat directories are scaffolded. |

## Troubleshooting

### Mapping looks wrong

- prefer MOL2 or SDF if PDB loses bond-order information
- verify protonation and tautomer states before building
- inspect `mapping.html` first

### Short run crashes during EM/NVT/NPT

- inspect `common/hybrid/hybrid.itp` and `common/hybrid/mapping.html`
- inspect `build/em.log`, `build/nvt.log`, and `build/npt.log`
- rerun a smoke test before attempting long production

### Analysis CLI finds no windows

- pass a repeat directory, not only the FEP root
- confirm that `window_*` directories exist and contain production outputs

## Summary

You have now:

- built a PRISM FEP scaffold for the 42-38 pair
- inspected the atom mapping and hybrid topology
- understood how PRISM currently transforms Coulomb and VDW terms across lambda windows
- run smoke-test and full-production entry points
- seen the main execution and analysis options exposed by the current code

## Additional Resources

- [FEP User Guide](../user-guide/fep-calculations.md)
- [Configuration](../user-guide/configuration.md)
- [Analysis Tools](../user-guide/analysis-tools.md)
- [Test Systems](https://github.com/AIB001/PRISM/tree/main/tests/gxf/FEP/unit_test)

## References

[^gmx-current]: GROMACS current manual. https://manual.gromacs.org/current/
[^gmx-free-energy-impl]: GROMACS current reference manual, *Free energy implementation*. https://manual.gromacs.org/current/reference-manual/special/free-energy-implementation.html
[^gmx-free-energy-interactions]: GROMACS current reference manual, *Free-energy interactions*. https://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html
[^gmx-grompp]: GROMACS current online help, `gmx grompp`. https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html
