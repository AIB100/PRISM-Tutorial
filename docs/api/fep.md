# FEP API

The PRISM FEP module provides comprehensive tools for Free Energy Perturbation calculations, including atom mapping, hybrid topology construction, system building, and multi-estimator analysis.

## Core Classes

### FEPScaffoldBuilder

Main class for building FEP systems with bound and unbound legs.

```python
from prism.fep import FEPScaffoldBuilder

builder = FEPScaffoldBuilder(
    protein_path="protein.pdb",
    ligand_a_path="ligand_a.mol2",
    ligand_b_path="ligand_b.mol2",
    output_dir="fep_output",
    forcefield="amber14sb",
    ligand_forcefield="gaff2",
    fep_config="fep.yaml"
)

# Build complete FEP scaffold
fep_dir = builder.build()
```

#### Parameters

- **protein_path** (str): Path to protein PDB file
- **ligand_a_path** (str): Path to reference ligand (state A)
- **ligand_b_path** (str): Path to mutant ligand (state B)
- **output_dir** (str): Output directory for FEP scaffold
- **forcefield** (str): Protein force field name
- **ligand_forcefield** (str): Ligand force field (gaff, gaff2, openff, cgenff, opls, rtf)
- **fep_config** (str): Path to FEP configuration YAML file
- **replicas** (int): Number of independent repeats (default: 1)

#### Methods

##### build()

Build complete FEP scaffold with bound and unbound legs.

```python
fep_dir = builder.build()
```

**Returns**: Path to FEP scaffold directory (`GMX_PROLIG_FEP/`)

**Output structure**:
```
GMX_PROLIG_FEP/
├── bound/
│   ├── repeat1/
│   │   ├── window_00/  # Lambda window directories
│   │   ├── window_01/
│   │   └── ...
│   └── repeat2/
├── unbound/
│   └── (same structure)
└── common/
    └── hybrid/
        ├── hybrid.itp
        ├── hybrid.gro
        └── mapping.html
```

##### build_bound_leg()

Build bound leg (protein-ligand complex) only.

```python
bound_dir = builder.build_bound_leg()
```

**Returns**: Path to bound leg directory

##### build_unbound_leg()

Build unbound leg (ligand in water) only.

```python
unbound_dir = builder.build_unbound_leg()
```

**Returns**: Path to unbound leg directory

---

### DistanceAtomMapper

Performs distance-based atom mapping between two ligands.

```python
from prism.fep import DistanceAtomMapper

mapper = DistanceAtomMapper(
    dist_cutoff=0.6,
    charge_cutoff=0.05,
    charge_strategy="mean"
)

# Perform mapping
mapping = mapper.map_atoms(
    ligand_a_atoms,
    ligand_b_atoms
)
```

#### Parameters

- **dist_cutoff** (float): Maximum distance for atom matching in Å (default: 0.6)
- **charge_cutoff** (float): Maximum charge difference for matching (default: 0.05)
- **charge_strategy** (str): Charge assignment for common atoms ("ref", "mut", "mean")
- **use_atom_types** (bool): Whether to consider atom types in matching (default: True)

#### Methods

##### map_atoms()

Map atoms between two ligands.

```python
mapping = mapper.map_atoms(
    atoms_a,  # List of Atom objects from ligand A
    atoms_b   # List of Atom objects from ligand B
)
```

**Returns**: `AtomMapping` object with classification results

**Atom classification**:
- **Common**: Present in both ligands with same position/charge
- **Transformed A**: Present only in ligand A (disappearing)
- **Transformed B**: Present only in ligand B (appearing)
- **Surrounding**: Near mutation site but not directly involved

##### get_mapping_statistics()

Get statistics about the mapping quality.

```python
stats = mapper.get_mapping_statistics()
print(stats)
# {'common': 28, 'transformed_a': 2, 'transformed_b': 2, 'surrounding': 0}
```

---

### HybridTopologyBuilder

Constructs single-topology hybrid ligand with A/B-state encoding.

```python
from prism.fep import HybridTopologyBuilder

hybrid_builder = HybridTopologyBuilder(
    mapping=atom_mapping,
    ligand_a_topology=topology_a,
    ligand_b_topology=topology_b
)

# Build hybrid topology
hybrid_topology = hybrid_builder.build()
```

#### Parameters

- **mapping** (AtomMapping): Atom mapping from DistanceAtomMapper
- **ligand_a_topology**: GROMACS topology for ligand A
- **ligand_b_topology**: GROMACS topology for ligand B
- **charge_reception** (str): How to handle charge redistribution ("surround", "none")

#### Methods

##### build()

Build hybrid topology with A/B-state encoding.

```python
hybrid_topology = hybrid_builder.build()
```

**Returns**: `HybridTopology` object with atoms and parameters

**Output format** (ITP file):
```
[ atoms ]
;   nr       type  resnr residue  atom   cgnr    charge       mass   typeB    chargeB    massB
   1        ca      1    LIG       CA     1     -0.120    12.011    ca      -0.120    12.011   ; Common
   2        ha      1    LIG       HA     1      0.150     1.008    ha       0.150     1.008   ; Common
   3        c       1    LIG       C1     1     -0.115    12.011    c       -0.115    12.011   ; A-only
   4        hc      1    LIG       H1     1      0.150     1.008   dummy      0.000    12.011   ; Transformed A
   5     dummy      1    LIG       C2     1      0.000    12.011     c        0.180    12.011   ; Transformed B
   6     dummy      1    LIG       H2     1      0.000     1.008    hc       0.150     1.008   ; B-only
```

---

## Analysis Classes

### FEPAnalyzer

Analyzes FEP simulation results with multiple estimators.

```python
from prism.fep.analysis import FEPAnalyzer

analyzer = FEPAnalyzer(
    bound_dir="fep_output/GMX_PROLIG_FEP/bound",
    unbound_dir="fep_output/GMX_PROLIG_FEP/unbound",
    temperature=310.0  # Kelvin
)

# Run analysis with MBAR estimator
results = analyzer.analyze(
    estimator="MBAR",
    n_bootstrap=1000,
    skip_frames=100  # Skip equilibration
)
```

#### Parameters

- **bound_dir** (str): Path to bound leg directory
- **unbound_dir** (str): Path to unbound leg directory
- **temperature** (float): Simulation temperature in Kelvin (default: 310.0)
- **repeats** (list): List of repeat numbers to analyze (default: all)

#### Methods

##### analyze()

Run FEP analysis with specified estimator.

```python
results = analyzer.analyze(
    estimator="MBAR",      # "BAR", "MBAR", or "TI"
    n_bootstrap=1000,      # Bootstrap iterations for uncertainty
    n_jobs=8,              # Parallel jobs for bootstrap
    skip_frames=100,       # Frames to skip for equilibration
    repeats=None           # Specific repeats (None = all)
)
```

**Returns**: Dictionary with analysis results

```python
{
    'bound_dg': -5.23,           # Bound leg ΔG (kcal/mol)
    'unbound_dg': -3.45,         # Unbound leg ΔG (kcal/mol)
    'binding_dg': 1.78,          # ΔΔG binding (kcal/mol)
    'bound_error': 0.15,         # Bound leg uncertainty
    'unbound_error': 0.12,       # Unbound leg uncertainty
    'binding_error': 0.20,       # ΔΔG uncertainty
    'overlap_matrix': [...],     # Overlap between windows
    'dhdl_profiles': [...]       # ∂H/∂λ profiles
}
```

##### generate_html_report()

Generate interactive HTML report.

```python
analyzer.generate_html_report(
    output_path="fep_results.html",
    include_plots=True,
    include_overlap=True,
    include_convergence=True
)
```

---

### FEPMultiEstimatorAnalyzer

Run multiple estimators and generate comparison report.

```python
from prism.fep.analysis import FEPMultiEstimatorAnalyzer

multi_analyzer = FEPMultiEstimatorAnalyzer(
    bound_dir="fep_output/GMX_PROLIG_FEP/bound",
    unbound_dir="fep_output/GMX_PROLIG_FEP/unbound",
    temperature=310.0
)

# Run all estimators
results = multi_analyzer.analyze_all(
    estimators=["BAR", "MBAR", "TI"],
    n_bootstrap=1000,
    n_jobs=8
)

# Generate comparison report
multi_analyzer.generate_html_report(
    output_path="fep_comparison.html"
)
```

#### Methods

##### analyze_all()

Run analysis with multiple estimators.

```python
results = multi_analyzer.analyze_all(
    estimators=["BAR", "MBAR", "TI"],
    n_bootstrap=1000,
    n_jobs=8,
    skip_frames=100
)
```

**Returns**: Dictionary with results from all estimators

```python
{
    'BAR': {'binding_dg': 1.75, 'binding_error': 0.18, ...},
    'MBAR': {'binding_dg': 1.78, 'binding_error': 0.20, ...},
    'TI': {'binding_dg': 1.82, 'binding_error': 0.22, ...}
}
```

##### compare_estimators()

Compare results across estimators.

```python
comparison = multi_analyzer.compare_estimators()
print(comparison)
# BAR:  1.75 ± 0.18 kcal/mol
# MBAR: 1.78 ± 0.20 kcal/mol
# TI:   1.82 ± 0.22 kcal/mol
#
# Consistency: All estimators agree within error bars
```

---

## Configuration

### FEPConfig

Load and manage FEP configuration from YAML files.

```python
from prism.fep import FEPConfig, read_fep_config

# Load from YAML file
config = read_fep_config("fep.yaml")

# Access configuration sections
print(config.mapping.dist_cutoff)        # 0.6
print(config.lambda.windows)            # 11
print(config.simulation.temperature)    # 310
```

#### Configuration Structure

```yaml
# Mapping parameters
mapping:
  dist_cutoff: 0.6
  charge_cutoff: 0.05
  charge_common: mean
  charge_reception: surround

# Lambda windows
lambda:
  strategy: decoupled     # or 'coupled'
  distribution: nonlinear  # or 'linear'
  windows: 11
  coul_windows: 4         # For decoupled strategy
  vdw_windows: 7          # For decoupled strategy

# Simulation parameters
simulation:
  temperature: 310        # Kelvin
  pressure: 1.0           # bar
  production_time_ns: 2.0 # Production per window
  dt: 0.002              # Timestep (ps)
  equilibration_nvt_time_ps: 500
  equilibration_npt_time_ps: 500

# Execution parameters
execution:
  mode: standard          # or 'repex'
  total_cpus: 56
  num_gpus: 4
  parallel_windows: 4
  use_gpu_pme: true
```

---

## Data Structures

### Atom

Represents an atom in a ligand.

```python
from prism.fep import Atom

atom = Atom(
    index=1,
    name="CA",
    type="ca",
    charge=-0.12,
    mass=12.011,
    coord=np.array([1.0, 2.0, 3.0])
)
```

#### Attributes

- **index** (int): Atom index (1-based)
- **name** (str): Atom name
- **type** (str): Atom type (force field specific)
- **charge** (float): Partial charge (e)
- **mass** (float): Atomic mass (amu)
- **coord** (np.ndarray): 3D coordinates (Å)

---

### HybridAtom

Represents an atom in hybrid topology with A/B states.

```python
from prism.fep import HybridAtom

hybrid_atom = HybridAtom(
    index=1,
    name="CA",
    type_a="ca",
    type_b="ca",
    charge_a=-0.12,
    charge_b=-0.12,
    mass_a=12.011,
    mass_b=12.011,
    classification="common",  # or 'transformed_a', 'transformed_b'
    coord=np.array([1.0, 2.0, 3.0])
)
```

#### Attributes

- **index** (int): Atom index
- **name** (str): Atom name
- **type_a** (str): Atom type in state A
- **type_b** (str): Atom type in state B
- **charge_a** (float): Charge in state A
- **charge_b** (float): Charge in state B
- **mass_a** (float): Mass in state A
- **mass_b** (float): Mass in state B
- **classification** (str): Atom classification
- **coord** (np.ndarray): 3D coordinates

---

### AtomMapping

Container for atom mapping results.

```python
mapping = mapper.map_atoms(atoms_a, atoms_b)

# Access mapping results
print(mapping.common_atoms)        # List of common atom pairs
print(mapping.transformed_a_atoms) # List of A-only atoms
print(mapping.transformed_b_atoms) # List of B-only atoms
print(mapping.surrounding_atoms)   # List of surrounding atoms
```

#### Methods

##### get_statistics()

Get mapping statistics.

```python
stats = mapping.get_statistics()
# {'common': 28, 'transformed_a': 2, 'transformed_b': 2, 'surrounding': 0}
```

##### validate()

Validate mapping quality.

```python
is_valid = mapping.validate()
# Returns True if:
# - No gray atoms (all classified)
# - Total charge ≈ 0 for neutral molecules
# - Reasonable number of transformed atoms
```

---

## CLI Interface

### FEP Analysis CLI

Command-line interface for analyzing FEP results.

```bash
python -m prism.fep.analysis.cli \
  --bound-dir fep_output/GMX_PROLIG_FEP/bound \
  --unbound-dir fep_output/GMX_PROLIG_FEP/unbound \
  --estimator BAR MBAR TI \
  --n-bootstrap 1000 \
  --n-jobs 8 \
  --output fep_results.html
```

#### Arguments

- `--bound-dir`: Path to bound leg directory (required)
- `--unbound-dir`: Path to unbound leg directory (required)
- `--estimator`: Estimator(s) to use (BAR, MBAR, TI)
- `--all-estimators`: Run all available estimators
- `--n-bootstrap`: Number of bootstrap iterations (default: 1000)
- `--n-jobs`: Parallel jobs for bootstrap (default: 4)
- `--skip-frames`: Frames to skip for equilibration (default: 0)
- `--repeats`: Specific repeats to analyze (default: all)
- `--output`: Output HTML file path (required)
- `--temperature`: Simulation temperature in Kelvin (default: 310)
- `--json`: Save results to JSON file (optional)

---

## Utility Functions

### Naming Functions

Generate standardized FEP system names.

```python
from prism.fep import generate_fep_system_name

system_name = generate_fep_system_name(
    protein_ff="amber14sb",
    ligand_ff="gaff2"
)
# Returns: "amber14sb-mut_gaff2"
```

### Validation Functions

Validate FEP system naming.

```python
from prism.fep import validate_fep_system_name

is_valid = validate_fep_system_name("amber14sb-mut_gaff2")
# Returns: True
```

---

## Complete Example

```python
from prism.fep import (
    FEPScaffoldBuilder,
    DistanceAtomMapper,
    FEPAnalyzer
)
from prism.fep import read_fep_config

# 1. Load configuration
config = read_fep_config("fep.yaml")

# 2. Build FEP system
builder = FEPScaffoldBuilder(
    protein_path="protein.pdb",
    ligand_a_path="ligand_a.mol2",
    ligand_b_path="ligand_b.mol2",
    output_dir="fep_output",
    forcefield="amber14sb",
    ligand_forcefield="gaff2",
    fep_config="fep.yaml",
    replicas=3
)

fep_dir = builder.build()
print(f"FEP scaffold built: {fep_dir}")

# 3. Run simulations (external)
# bash fep_output/GMX_PROLIG_FEP/run_fep.sh all

# 4. Analyze results
analyzer = FEPAnalyzer(
    bound_dir=f"{fep_dir}/bound",
    unbound_dir=f"{fep_dir}/unbound",
    temperature=310.0
)

results = analyzer.analyze(
    estimator="MBAR",
    n_bootstrap=1000,
    n_jobs=8
)

print(f"ΔΔG binding: {results['binding_dg']:.2f} ± {results['binding_error']:.2f} kcal/mol")

# 5. Generate report
analyzer.generate_html_report("fep_results.html")
```

---

## Error Handling

### Common Exceptions

```python
from prism.exceptions import (
    FEPMappingError,
    FEPTopologyError,
    FEPAnalysisError
)

try:
    mapping = mapper.map_atoms(atoms_a, atoms_b)
except FEPMappingError as e:
    print(f"Mapping failed: {e}")
    # Try adjusting cutoffs
    mapper.dist_cutoff = 0.8
    mapping = mapper.map_atoms(atoms_a, atoms_b)

try:
    hybrid = hybrid_builder.build()
except FEPTopologyError as e:
    print(f"Topology construction failed: {e}")
    # Check for missing B-state parameters

try:
    results = analyzer.analyze(estimator="MBAR")
except FEPAnalysisError as e:
    print(f"Analysis failed: {e}")
    # Check for missing dhdl.xvg files
```

---

## Performance Tips

1. **Parallelize lambda windows**: Use `parallel_windows` in config
2. **Bootstrap parallelization**: Set `n_jobs` in analysis
3. **Memory-efficient analysis**: Process one repeat at a time
4. **GPU acceleration**: Enable GPU for PME calculations
5. **Optimal window count**: Balance accuracy vs computational cost

---

## See Also

- [FEP User Guide](../user-guide/fep-calculations.md) - Complete FEP documentation
- [FEP Tutorial](../tutorials/fep-tutorial.md) - Step-by-step workflow
- [Builder API](builder.md) - General system building
- [Analysis API](analysis.md) - Trajectory analysis tools
