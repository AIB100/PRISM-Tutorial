# API Reference

Comprehensive API documentation for all PRISM modules.

## Core Modules

### System Building
- [PRISMBuilder](builder.md) - Main system builder class
- [PRISMSystem](builder.md#prismsystem) - High-level system interface

### Analysis
- [TrajAnalysis](analysis.md) - Trajectory analysis tools
- [Visualization](analysis.md#visualization) - Interactive visualizations

### Advanced Calculations
- [PMF Module](pmf.md) - Free energy calculations
- [Utilities](utilities.md) - Helper functions

---

## Quick Reference

### High-Level API

```python
import prism

# Build system
system = prism.system("protein.pdb", "ligand.mol2")
system.build()

# Analyze trajectory
analysis = prism.analyze_trajectory("system.gro", "traj.xtc")

# PMF calculation
results = prism.run_pmf_workflow("md_dir", "pmf_output")
```

### Class-Based API

```python
from prism import PRISMBuilder, TrajAnalysis

# Builder
builder = PRISMBuilder("protein.pdb", "ligand.mol2", "output")
builder.run()

# Analysis
analyzer = TrajAnalysis("topology.gro", "trajectory.xtc")
analyzer.analyze_all()
```

---

## Module Index

| Module | Purpose | Key Classes |
|--------|---------|-------------|
| [builder](builder.md) | System construction | PRISMBuilder, PRISMSystem |
| [analysis](analysis.md) | Trajectory analysis | TrajAnalysis, HTMLGenerator |
| [pmf](pmf.md) | Free energy | PMFBuilder, PMFRunner |
| [utilities](utilities.md) | Helper functions | Various utility functions |

---

## API Stability

- **Stable**: High-level functions (prism.system, prism.analyze_trajectory)
- **Evolving**: Internal classes (subject to change)
- **Experimental**: PMF module (API may change in minor versions)
