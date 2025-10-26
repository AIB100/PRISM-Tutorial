# PMF Calculation Tutorial

## Overview

Calculate protein-ligand binding free energies using umbrella sampling and WHAM analysis.

**Time:** 2-4 hours + computation time
**Level:** Advanced

## Prerequisites

- Equilibrated MD system from [Basic Tutorial](basic-tutorial.md)
- Understanding of free energy concepts
- Access to computational resources (preferably GPUs)

## Workflow

1. Equilibrate protein-ligand complex
2. Steered MD to pull ligand away
3. Umbrella sampling along unbinding pathway
4. WHAM analysis to calculate PMF
5. Extract binding free energy

## Step 1: Prepare Equilibrated System

```bash
# Start from basic tutorial output
cd t4_lysozyme_benzene
```

## Step 2: One-Step PMF Calculation

```python
import prism as pm

# Automated PMF workflow
results = pm.run_pmf_workflow(
    md_system_dir=".",
    output_dir="pmf_results"
)

print(f"Binding energy: {results['binding_energy']['value']:.2f} ± {results['binding_energy']['error']:.2f} kcal/mol")
```

## Step 3: Step-by-Step PMF (Advanced)

### Generate SMD

```python
import prism as pm

# Create PMF system
pmf_sys = pm.pmf_system(".", "pmf_manual")

# Build SMD
smd_result = pmf_sys.build(step='smd')
print("SMD prepared in: pmf_manual/pmf_smd/")
```

Run SMD:
```bash
cd pmf_manual/pmf_smd
bash run_smd.sh
cd ../..
```

### Generate Umbrella Windows

```python
# Build umbrella sampling windows
umbrella_result = pmf_sys.build_umbrella_step()
print(f"Generated {umbrella_result['n_windows']} windows")
```

Run umbrella sampling:
```bash
cd pmf_manual/pmf_umbrella
bash run_all_umbrella.sh  # Or submit to cluster
cd ../..
```

### WHAM Analysis

```python
# Analyze with WHAM
analysis_result = pmf_sys.build_analysis_step()

print(f"ΔG_bind = {analysis_result['binding_energy']['value']:.2f} kcal/mol")
print(f"Error = {analysis_result['binding_energy']['error']:.2f} kcal/mol")
```

## Step 4: Visualize PMF Profile

```python
import matplotlib.pyplot as plt
import numpy as np

# Load PMF data
data = np.loadtxt("pmf_results/pmf_analysis/pmf_profile.dat")
distance = data[:, 0]
pmf = data[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(distance, pmf, 'b-', linewidth=2)
plt.xlabel('Distance (nm)', fontsize=14)
plt.ylabel('PMF (kJ/mol)', fontsize=14)
plt.title('Protein-Ligand Unbinding PMF')
plt.grid(True, alpha=0.3)
plt.axhline(y=0, color='k', linestyle='--', alpha=0.5)
plt.savefig('pmf_profile.png', dpi=300)
```

## Expected Results

For T4 Lysozyme-Benzene:
- **Binding energy**: -4 to -6 kcal/mol (experimental: ~-5 kcal/mol)
- **Binding site depth**: ~2-3 nm
- **PMF should show clear minimum** at binding site

## Troubleshooting

### Ligand Escapes Too Quickly

Reduce SMD pull rate:
```python
config = {'smd': {'pull_rate': 0.005}}
results = pm.run_pmf_workflow(".", "pmf_slow", config=config)
```

### Poor Window Overlap

Increase number of windows:
```python
config = {'umbrella': {'n_windows': 50}}
```

### Large Errors

Extend sampling time:
```python
config = {'umbrella': {'production_time_ps': 20000}}
```

## Next Steps

- Try different pulling directions
- Calculate PMF for multiple force fields
- Compare with experimental binding affinity

## References

See [PMF Calculations Guide](../user-guide/pmf-calculations.md) for detailed information.

