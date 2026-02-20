# Custom Configuration Example

!!! example "Quick Start"
    ```bash
    prism --export-config my_config.yaml
    prism protein.pdb ligand.mol2 -o output --config my_config.yaml
    ```

Advanced system setup using YAML configuration files for fine-grained control over all PRISM parameters.

## Overview

Learn how to:
- Create and customize YAML configuration files
- Control all aspects of system building
- Set advanced MD parameters
- Use configuration templates

**Complexity**: Advanced  
**Prerequisites**: Completed [Simple Example](simple.md)

---

## Configuration File Structure

### Complete Configuration Template

```yaml
# PRISM Configuration File
# All parameters with defaults shown

general:
  overwrite: false
  verbose: true

# Force fields
forcefields:
  protein: amber99sb
  ligand: gaff
  water: tip3p

# Box configuration
box:
  distance: 1.5  # nm from protein to box edge
  shape: cubic   # cubic, dodecahedron, octahedron
  center: true

# Solvation and ions
ions:
  neutral: true
  concentration: 0.15  # M
  positive_ion: NA
  negative_ion: CL

# Simulation parameters
simulation:
  temperature: 310  # K
  pressure: 1.0     # bar
  pH: 7.4
  ligand_charge: 0
  production_time_ns: 500
  dt: 0.002  # ps
  equilibration_nvt_time_ps: 500
  equilibration_npt_time_ps: 500

# Protein preparation
protein_preparation:
  ion_handling_mode: smart  # smart, keep_all, remove_all
  metal_distance_cutoff: 5.0
  keep_crystal_water: false
  remove_crystallization_artifacts: true

# Output options
output:
  trajectory_interval_ps: 500
  energy_interval_ps: 10
  compressed_trajectory: true

# Constraints
constraints:
  algorithm: lincs
  type: h-bonds
```

---

## Example Configurations

### 1. High-Accuracy Configuration

```yaml
# config_accurate.yaml
forcefields:
  protein: amber14sb
  ligand: openff
  water: tip4p

box:
  distance: 2.0
  shape: dodecahedron

simulation:
  temperature: 300
  production_time_ns: 1000
  dt: 0.001  # Smaller timestep
  equilibration_nvt_time_ps: 1000
  equilibration_npt_time_ps: 1000

output:
  trajectory_interval_ps: 100  # More frequent output
```

### 2. Fast Screening Configuration

```yaml
# config_fast.yaml
box:
  distance: 1.0  # Minimal box
  shape: cubic

simulation:
  production_time_ns: 50
  equilibration_nvt_time_ps: 250
  equilibration_npt_time_ps: 250

output:
  trajectory_interval_ps: 1000
  compressed_trajectory: true
```

### 3. Membrane Protein Configuration

```yaml
# config_membrane.yaml
forcefields:
  protein: charmm36
  ligand: cgenff
  water: tips3p

simulation:
  temperature: 310
  pressure: 1.0

# No box for membrane (use pre-built membrane system)
```

---

## Usage

```python
import prism

# Use custom configuration
system = prism.system(
    "protein.pdb",
    "ligand.mol2",
    output_dir="custom_system",
    config="config_accurate.yaml"
)

system.build()
```

---

## Advanced Features

### Dynamic Parameter Overrides

```python
# Start with config file, override specific parameters
system = prism.system(
    "protein.pdb",
    "ligand.mol2",
    config="config_accurate.yaml",
    production_ns=2000,  # Override config value
    temperature=320      # Override config value
)
```

### Configuration Validation

```python
from prism.utils.config import ConfigurationManager

config = ConfigurationManager("my_config.yaml")
config.validate()  # Check for errors
config.show()      # Display all settings
```

---

## Best Practices

1. **Start with Template**: Export default config and modify
2. **Version Control**: Keep configs in git
3. **Document Changes**: Add comments to explain non-standard settings
4. **Validate**: Test configs on small systems first

---

<div class="whats-next" markdown>

## What's Next

- [Read the full Configuration Guide](../user-guide/configuration.md)
- [See the Builder API Reference](../api/builder.md)
- [Return to Examples overview](index.md)

</div>
