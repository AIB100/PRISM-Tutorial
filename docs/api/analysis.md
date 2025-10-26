# Analysis API

Trajectory analysis and visualization tools.

## TrajAnalysis

Main class for MD trajectory analysis.

```python
from prism import TrajAnalysis

analyzer = TrajAnalysis(
    topology="system.gro",
    trajectory="md.xtc",
    ligand_resname="LIG"
)
```

### Parameters

- **topology** (str): Topology file (.gro, .pdb, .tpr)
- **trajectory** (str): Trajectory file (.xtc, .dcd, .trr)
- **ligand_resname** (str): Ligand residue name (default: "LIG")

### Methods

#### analyze_all()
Run complete analysis suite.

```python
analyzer.analyze_all(output_dir="analysis_results")
```

**Generates**:
- RMSD plots
- RMSF per-residue
- Contact maps
- H-bond analysis
- Interactive HTML

#### calc_rmsd()
Calculate RMSD.

```python
rmsd = analyzer.calc_rmsd(
    selection="protein",
    reference_frame=0
)
```

**Returns**: numpy array of RMSD values (Å)

#### calc_rmsf()
Calculate RMSF.

```python
rmsf = analyzer.calc_rmsf(selection="protein")
```

**Returns**: numpy array of RMSF values (Å)

#### calc_contacts()
Calculate protein-ligand contacts.

```python
contacts = analyzer.calc_contacts(
    cutoff=3.5,  # Å
    scheme='closest-heavy'
)
```

**Returns**: Contact frequency array

#### calc_hbonds()
Calculate hydrogen bonds.

```python
hbonds = analyzer.calc_hbonds(
    distance_cutoff=3.5,  # Å
    angle_cutoff=120  # degrees
)
```

**Returns**: H-bond data

#### calc_sasa()
Calculate solvent accessible surface area.

```python
sasa = analyzer.calc_sasa(selection="resname LIG")
```

**Returns**: SASA values (nm²)

---

## Visualization

### visualize_trajectory()

Generate interactive HTML visualization.

```python
import prism

html_file = prism.visualize_trajectory(
    trajectory="md.xtc",
    topology="system.gro",
    ligand="ligand.sdf",
    output="contacts.html"
)
```

**Parameters**:
- trajectory (str): Trajectory file
- topology (str): Topology file
- ligand (str): Ligand structure file
- output (str): Output HTML file

**Returns**: Path to HTML file

### HTMLGenerator

Advanced visualization control.

```python
from prism.analysis.visualization import HTMLGenerator

generator = HTMLGenerator(
    "trajectory.xtc",
    "topology.gro",
    "ligand.sdf"
)

generator.analyze()
generator.generate("output.html")
```

---

## Trajectory Processing

### process_trajectory()

Fix PBC artifacts and center system.

```python
import prism

processed = prism.process_trajectory(
    input_trajectory="raw.dcd",
    output_trajectory="processed.xtc",
    topology_file="system.tpr",
    center_selection="Protein",
    pbc_method="mol"
)
```

### batch_process_trajectories()

Process multiple trajectories.

```python
trajectories = ["traj1.dcd", "traj2.dcd", "traj3.dcd"]

processed = prism.batch_process_trajectories(
    input_trajectories=trajectories,
    output_dir="processed",
    topology_file="system.tpr"
)
```

---

## Examples

### Basic Analysis

```python
from prism import TrajAnalysis

analyzer = TrajAnalysis("system.gro", "md.xtc")

# RMSD
rmsd_protein = analyzer.calc_rmsd(selection="protein")
rmsd_ligand = analyzer.calc_rmsd(selection="resname LIG")

# Contacts
contacts = analyzer.calc_contacts()
print(f"Average contacts: {contacts.mean()}")
```

### Complete Workflow

```python
import prism

# Analyze trajectory
analyzer = prism.analyze_trajectory(
    "system.gro",
    "md.xtc",
    ligand_resname="LIG",
    output_dir="analysis"
)

# Generate visualization
prism.visualize_trajectory(
    "md.xtc",
    "system.gro",
    "ligand.sdf",
    output="interactive_contacts.html"
)
```

---

## Related

- [User Guide: Analysis Tools](../user-guide/analysis-tools.md)
- [Examples](../examples/simple.md)
