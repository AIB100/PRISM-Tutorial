# FEP API

The PRISM FEP module provides tools for:

- distance-based atom mapping
- hybrid topology construction with A/B-state parameters
- scaffold generation for bound and unbound legs
- FEP post-processing for `window_*` outputs

## Core Classes

### FEPScaffoldBuilder

`FEPScaffoldBuilder` stages a PRISM-style FEP workspace around an existing hybrid ligand package.

The most direct entry point is `build(receptor_pdb, hybrid_ligand_dir)`. For PRISM-generated ligand force-field components, use `build_from_components(...)`.

```python
from prism.fep.modeling import FEPScaffoldBuilder

builder = FEPScaffoldBuilder(
    output_dir="fep_output",
    config={"replicas": 1},
    lambda_windows=32,
    lambda_strategy="decoupled",
    lambda_distribution="nonlinear",
    overwrite=False,
)

fep_dir = builder.build(
    receptor_pdb="protein.pdb",
    hybrid_ligand_dir="hybrid_ligand_dir",
)
```

#### Constructor Parameters

- `output_dir` (`str`): output scaffold directory
- `config` (`dict | None`): scaffold/build configuration
- `lambda_windows` (`int`): number of lambda windows
- `lambda_strategy` (`str`): typically `decoupled` or `coupled`
- `lambda_distribution` (`str`): lambda spacing strategy
- `overwrite` (`bool`): whether to overwrite existing files

#### Methods

##### `build(...)`

Build the complete FEP scaffold.

**Returns**: `Path` pointing to `GMX_PROLIG_FEP/`

Typical layout:

```text
GMX_PROLIG_FEP/
├── run_fep.sh
├── fep_scaffold.json
├── common/
│   └── hybrid/
│       ├── hybrid.itp
│       ├── hybrid.gro
│       └── mapping.html
├── bound/
│   └── repeat1/
│       ├── build/
│       ├── input/
│       ├── mdps/
│       ├── run_prod_standard.sh
│       ├── run_prod_repex.sh
│       └── window_00/ ... window_N/
└── unbound/
    └── repeat1/
        └── window_00/ ... window_N/
```

### DistanceAtomMapper

`DistanceAtomMapper` performs distance-based atom mapping between ligand A and ligand B.

```python
from prism.fep import DistanceAtomMapper

mapper = DistanceAtomMapper(
    dist_cutoff=0.6,
    charge_cutoff=0.05,
    charge_common="mean",
    charge_reception="surround",
)

mapping = mapper.map(ligand_a_atoms, ligand_b_atoms)
```

#### Parameters

- `dist_cutoff` (`float`): atom-matching cutoff in **nm**
- `charge_cutoff` (`float`): charge-difference threshold
- `charge_common` (`str`): `ref`, `mut`, `mean`, or `none`
- `charge_reception` (`str`): `unique`, `surround`, or `none`
- `recharge_hydrogen` (`bool`): whether hydrogen charges can be perturbed

#### Notes

- Matching candidates require the same element and a distance below `dist_cutoff`.
- Generic/sequential atom types such as OpenFF, OPLS/LigParGen, and SwissParam relax type-based compatibility checks.
- GAFF/GAFF2 remains type-aware.

#### Methods

##### `map(ligand_a, ligand_b)`

Returns an `AtomMapping` object containing the final classification.

##### `from_config(config)`

Constructs a mapper from a loaded PRISM config dictionary.

### HybridTopologyBuilder

Builds the ligand hybrid topology using the atom mapping and the A/B ligand topologies.

```python
from prism.fep import HybridTopologyBuilder

hybrid_builder = HybridTopologyBuilder(
    mapping=atom_mapping,
    ligand_a_topology=topology_a,
    ligand_b_topology=topology_b,
    charge_reception="surround",
)

hybrid_topology = hybrid_builder.build()
```

The generated hybrid topology stores A-state and B-state parameters in the same molecular description.[^gmx-free-energy-impl] [^gmx-free-energy-interactions]

## Analysis Classes

### FEPAnalyzer

Analyzes one set of bound and unbound FEP windows.

```python
from prism.fep.analysis import FEPAnalyzer

analyzer = FEPAnalyzer(
    bound_dir="fep_output/GMX_PROLIG_FEP/bound/repeat1",
    unbound_dir="fep_output/GMX_PROLIG_FEP/unbound/repeat1",
    temperature=310.0,
    estimator="MBAR",
    bootstrap_n_jobs=4,
)

results = analyzer.analyze()
report_path = analyzer.generate_html_report("fep_results.html")
```

#### Parameters

- `bound_dir` (`str | Path | list`): bound repeat directory or repeat directories
- `unbound_dir` (`str | Path | list`): unbound repeat directory or repeat directories
- `temperature` (`float`): simulation temperature in Kelvin
- `estimator` (`str`): `MBAR`, `BAR`, or `TI`
- `backend` (`str`): `alchemlyb` or `gmx_bar`
- `energy_components` (`list[str] | None`): defaults to `['elec', 'vdw']`
- `bootstrap_n_jobs` (`int`): bootstrap worker count

##### `analyze()`

Runs the configured estimator and returns a `FEResults` object.

##### `generate_html_report(output_path)`

Generates an HTML report for the current results.

### FEPMultiEstimatorAnalyzer

Runs multiple estimators and compares them in one report.

```python
from prism.fep.analysis import FEPMultiEstimatorAnalyzer

multi = FEPMultiEstimatorAnalyzer(
    bound_dirs="fep_output/GMX_PROLIG_FEP/bound/repeat1",
    unbound_dirs="fep_output/GMX_PROLIG_FEP/unbound/repeat1",
    estimators=["TI", "BAR", "MBAR"],
    temperature=310.0,
    bootstrap_n_jobs=4,
)

results = multi.analyze()
```

#### Parameters

- `bound_dirs` (`str | Path | list`): bound repeat directory or directories
- `unbound_dirs` (`str | Path | list`): unbound repeat directory or directories
- `estimators` (`list[str] | None`): estimator list
- `temperature` (`float`): simulation temperature in Kelvin
- `backend` (`str`): currently `alchemlyb` for multi-estimator mode
- `bootstrap_n_jobs` (`int`): bootstrap worker count

## Configuration

### FEPConfig

`FEPConfig` loads and merges PRISM FEP configuration from YAML files.

```python
from prism.fep.common.config import FEPConfig

cfg = FEPConfig(work_dir=".")
print(cfg.get_mapping_params())
print(cfg.get_lambda_params())
print(cfg.get_simulation_params())
```

Current canonical YAML structure:

```yaml
mapping:
  dist_cutoff: 0.6        # nm
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
  temperature: 310
  pressure: 1.0
  production_time_ns: 2.0
  dt: 0.002
  equilibration_nvt_time_ps: 500
  equilibration_npt_time_ps: 500

execution:
  mode: standard
  total_cpus: 56
  num_gpus: 4
  parallel_windows: 4
  use_gpu_pme: true
```

### `read_fep_config(config_file)`

Reads the legacy FEbuilder-compatible config format and returns a plain dictionary.

## Data Structures

### Atom

Represents a ligand atom.

Key fields:

- `index`
- `name`
- `type`
- `charge`
- `mass`
- `coord`

Coordinates used by the current mapper are interpreted in **nm**.

### HybridAtom

Represents a hybrid-topology atom with A/B-state parameters.

Key fields:

- `type_a`, `type_b`
- `charge_a`, `charge_b`
- `mass_a`, `mass_b`
- `classification`

### AtomMapping

Container for the final atom-mapping classification.

Typical access patterns:

```python
mapping = mapper.map(atoms_a, atoms_b)
print(mapping.common_atoms)
print(mapping.transformed_a_atoms)
print(mapping.transformed_b_atoms)
print(mapping.surrounding_atoms)
```

## CLI Interface

### FEP Analysis CLI

Current supported CLI entry points are:

```bash
prism --fep-analyze \
  --bound-dir fep_output/GMX_PROLIG_FEP/bound/repeat1 \
  --unbound-dir fep_output/GMX_PROLIG_FEP/unbound/repeat1 \
  --estimator MBAR BAR TI \
  --bootstrap-n-jobs 8 \
  --output fep_results.html \
  --json fep_results.json
```

Equivalent module entry point:

```bash
python -m prism.fep.analysis.cli \
  --bound-dir fep_output/GMX_PROLIG_FEP/bound/repeat1 \
  --unbound-dir fep_output/GMX_PROLIG_FEP/unbound/repeat1 \
  --estimator MBAR BAR TI \
  --bootstrap-n-jobs 8 \
  --output fep_results.html
```

#### Arguments

- `--bound-dir`: repeat directory containing `window_*`
- `--unbound-dir`: repeat directory containing `window_*`
- `--estimator`: one or more of `BAR`, `MBAR`, `TI`
- `--all-estimators`: run all estimators
- `--bootstrap-n-jobs`: bootstrap worker count
- `--output`: HTML report path
- `--json`: optional JSON output
- `--temperature`: simulation temperature in Kelvin
- `--backend`: analysis backend

## Utility Functions

### Naming Functions

```python
from prism.fep import generate_fep_system_name, validate_fep_system_name

name = generate_fep_system_name("amber14sb_OL15", "gaff2")
assert validate_fep_system_name(name)
```

## References

[^gmx-current]: GROMACS current manual. https://manual.gromacs.org/current/
[^gmx-free-energy-impl]: GROMACS current reference manual, *Free energy implementation*. https://manual.gromacs.org/current/reference-manual/special/free-energy-implementation.html
[^gmx-free-energy-interactions]: GROMACS current reference manual, *Free-energy interactions*. https://manual.gromacs.org/current/reference-manual/functions/free-energy-interactions.html
[^gmx-grompp]: GROMACS current online help, `gmx grompp`. https://manual.gromacs.org/current/onlinehelp/gmx-grompp.html
[^fep-user-guide]: PRISM Tutorial, *FEP Calculations*. ../user-guide/fep-calculations.md
[^fep-tutorial]: PRISM Tutorial, *FEP Workflow Tutorial*. ../tutorials/fep-tutorial.md
