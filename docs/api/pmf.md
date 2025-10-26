# PMF API

Free energy calculations using umbrella sampling and WHAM.

## PMFBuilder

System preparation for PMF calculations.

```python
from prism.pmf import PMFBuilder

builder = PMFBuilder(
    md_results_dir="md_system/GMX_PROLIG_MD",
    output_dir="pmf_system"
)

builder.build(equilibrate=True)
```

### Methods

#### build()
Build PMF-optimized system.

```python
builder.build(
    equilibrate=True,
    extend_box=True
)
```

**Parameters**:
- equilibrate (bool): Run equilibration
- extend_box (bool): Extend box for PMF

---

## PMFRunner

Complete PMF workflow execution.

```python
from prism.pmf import PMFRunner

runner = PMFRunner(config="pmf_config.yaml")

results = runner.run_complete_workflow(
    system_dir="pmf_system",
    output_dir="pmf_results"
)
```

### Methods

#### run_smd()
Run steered molecular dynamics.

```python
smd_results = runner.run_smd(
    system_dir="pmf_system",
    output_dir="pmf_smd",
    pull_rate=0.01,
    pull_distance=4.0
)
```

#### setup_umbrella_sampling()
Setup umbrella windows.

```python
umbrella_results = runner.setup_umbrella_sampling(
    smd_dir="pmf_smd",
    output_dir="umbrella",
    window_spacing=0.1,
    num_windows=40
)
```

#### analyze_pmf()
Perform WHAM analysis.

```python
wham_results = runner.analyze_pmf(
    umbrella_dir="umbrella",
    output_dir="analysis",
    temperature=310,
    num_bootstrap=100
)
```

---

## High-Level Functions

### run_pmf_workflow()

Complete automated PMF calculation.

```python
import prism

results = prism.run_pmf_workflow(
    md_system_dir="md_results",
    output_dir="pmf_output",
    config="pmf_config.yaml"
)

# Access results
binding_energy = results['binding_energy']['value']
error = results['binding_energy']['error']
```

**Returns**: Dictionary with:
- binding_energy: dict with 'value' and 'error'
- pmf_profile: PMF curve data
- analysis_files: paths to output files

### create_pmf_config()

Generate configuration template.

```python
import prism

prism.create_pmf_config(
    "my_pmf_config.yaml",
    template="accurate"  # or "fast", "default"
)
```

---

## Examples

### Automated PMF

```python
import prism

results = prism.run_pmf_workflow(
    "md_system",
    "pmf_results"
)

print(f"ΔG = {results['binding_energy']['value']:.2f} kcal/mol")
```

### Step-by-Step PMF

```python
from prism.pmf import PMFBuilder, PMFRunner

# 1. Prepare system
builder = PMFBuilder("md_system/GMX_PROLIG_MD", "pmf_system")
builder.build(equilibrate=True)

# 2. Run SMD
runner = PMFRunner()
smd = runner.run_smd("pmf_system", "smd_output")

# 3. Setup umbrella
umbrella = runner.setup_umbrella_sampling("smd_output", "umbrella")

# 4. Run umbrella simulations (manual)
# bash umbrella/run_all_windows.sh

# 5. Analyze
results = runner.analyze_pmf("umbrella", "analysis")
```

---

## Configuration

PMF configuration file structure:

```yaml
smd:
  pull_rate: 0.01  # nm/ps
  pull_distance: 4.0  # nm
  nsteps: 2000000
  force_constant: 1000  # kJ/mol/nm²

umbrella:
  window_spacing: 0.1  # nm
  num_windows: 40
  production_time_ps: 15000
  force_constant: 3000  # kJ/mol/nm²

wham:
  tolerance: 1e-6
  num_bootstrap: 100
  temperature: 310  # K
```

---

## Related

- [Tutorial: PMF Calculations](../tutorials/pmf-tutorial.md)
- [User Guide: PMF](../user-guide/pmf-calculations.md)
