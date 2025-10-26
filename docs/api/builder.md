# Builder API

The PRISMBuilder module provides automated protein-ligand system construction for GROMACS simulations.

## PRISMBuilder

Main class for building protein-ligand MD systems.

```python
from prism import PRISMBuilder

builder = PRISMBuilder(
    protein_path="protein.pdb",
    ligand_path="ligand.mol2",
    output_dir="my_system",
    ligand_forcefield="gaff",
    config_path=None,
    forcefield="amber99sb",
    water_model="tip3p"
)
```

### Parameters

- **protein_path** (str): Path to protein PDB file
- **ligand_path** (str): Path to ligand MOL2/SDF file
- **output_dir** (str): Output directory
- **ligand_forcefield** (str): Ligand force field (gaff, gaff2, openff, opls, etc.)
- **config_path** (str, optional): YAML configuration file
- **forcefield** (str, optional): Protein force field
- **water_model** (str, optional): Water model (tip3p, tip4p, spce)

### Methods

#### run()
Build complete MD system.

```python
builder.run()
```

#### generate_ligand_forcefield()
Generate ligand force field parameters.

```python
lig_ff_dir = builder.generate_ligand_forcefield()
```

**Returns**: Path to ligand force field directory

#### clean_protein()
Clean and prepare protein structure.

```python
cleaned = builder.clean_protein(
    ion_mode='smart',
    distance_cutoff=5.0,
    keep_crystal_water=False
)
```

**Parameters**:
- ion_mode (str): 'smart', 'keep_all', or 'remove_all'
- distance_cutoff (float): Distance in Å for metal ion filtering
- keep_crystal_water (bool): Keep crystallographic water

**Returns**: Path to cleaned protein file

#### build_model()
Build GROMACS topology and coordinate files.

```python
model_dir = builder.build_model(protein_file)
```

#### generate_mdp_files()
Generate MD parameter files.

```python
builder.generate_mdp_files()
```

#### generate_localrun_script()
Create automation script for running MD.

```python
script_path = builder.generate_localrun_script()
```

---

## PRISMSystem

High-level system interface.

```python
import prism

system = prism.system(
    "protein.pdb",
    "ligand.mol2",
    output_dir="output",
    ligand_forcefield="gaff",
    **kwargs
)
```

### Methods

#### build()
Build complete system.

```python
output_path = system.build()
```

**Returns**: Path to output directory

#### info()
Display system information.

```python
system.info()
```

#### get_output_files()
Get dictionary of all output file paths.

```python
files = system.get_output_files()
print(files['system_gro'])
print(files['system_top'])
```

**Returns**: dict with file paths

---

## Examples

### Basic Usage

```python
from prism import PRISMBuilder

builder = PRISMBuilder(
    "protein.pdb",
    "ligand.mol2",
    "output"
)
builder.run()
```

### Step-by-Step

```python
builder = PRISMBuilder("protein.pdb", "ligand.mol2", "output")

# Step 1: Ligand parameters
lig_ff = builder.generate_ligand_forcefield()

# Step 2: Clean protein
clean_prot = builder.clean_protein()

# Step 3: Build model
model = builder.build_model(clean_prot)

# Step 4: Generate MDPs
builder.generate_mdp_files()

# Step 5: Create run script
builder.generate_localrun_script()
```

### Custom Configuration

```python
builder = PRISMBuilder(
    "protein.pdb",
    "ligand.mol2",
    "output",
    config_path="custom_config.yaml"
)
builder.run()
```

---

## Related

- [Examples](../examples/simple.md)
- [User Guide: Building Systems](../user-guide/building-systems.md)
