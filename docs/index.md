---
hide:
  - toc
---

# PRISM

<div style="text-align: center;" markdown>

![PRISM Logo](assets/logo.png){ width="280" }

**Protein Receptor Interaction Simulation Modeler**

*A comprehensive toolkit for building and analyzing protein-ligand MD systems*

[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Docs](https://img.shields.io/badge/Docs-Latest-brightgreen.svg)](https://prism-docs.github.io)
[![GitHub](https://img.shields.io/badge/GitHub-PRISM-black.svg)](https://github.com/AIB001/PRISM)

</div>

---

## Overview

**PRISM** is a Python framework for **building protein-ligand systems** and **analyzing MD trajectories**, developed at the Institute of Quantitative Biology, Zhejiang University. It automates the complex workflow of preparing molecular dynamics simulations with multiple force field options and provides comprehensive analysis tools.

## Quick Start

=== "Command Line"

    ```bash
    # GAFF (default)
    prism protein.pdb ligand.mol2 -o my_system

    # OpenFF
    prism protein.pdb ligand.sdf -o my_system --ligand-forcefield openff

    # Multiple ligands
    prism -pf protein.pdb -lf lig1.mol2 -lf lig2.mol2 -o my_system
    ```

=== "Python API"

    ```python
    import prism

    # Build system with GAFF force field (default)
    system = prism.system("protein.pdb", "ligand.mol2", output_dir="my_system")
    system.build()
    ```

=== "Installation"

    ```bash
    git clone https://github.com/AIB001/PRISM.git
    cd PRISM
    pip install -e .
    ```

## Key Features

<div class="feature-grid" markdown>

<div class="feature-card" markdown>

### Automated System Building
One-command setup for protein-ligand MD systems with automated protein cleaning, ligand parameterization, solvation, and neutralization.

</div>

<div class="feature-card" markdown>

### Multiple Force Fields
Support for 8+ ligand force fields: GAFF, GAFF2, OpenFF, CGenFF, OPLS-AA, MMFF, MATCH, and hybrid approaches.

</div>

<div class="feature-card" markdown>

### Trajectory Analysis
Comprehensive analysis tools including RMSD/RMSF, contacts, hydrogen bonds, clustering, SASA, and dihedral angles.

</div>

<div class="feature-card" markdown>

### Interactive Visualization
Generate interactive HTML visualizations of protein-ligand contacts with customizable 2D/3D views.

</div>

<div class="feature-card" markdown>

### GROMACS Integration
Native GROMACS workflow with auto-generated MDP files for EM, NVT, NPT, and production MD.

</div>

<div class="feature-card" markdown>

### Advanced Calculations
Built-in support for PMF (umbrella sampling/WHAM) and MM/PBSA binding energy calculations.

</div>

</div>

## Core Modules

| Module | Description |
| --- | --- |
| **Builder** | Automated system construction: force field generation, protein preparation, solvation, topology output |
| **Simulation** | Complete EM → NVT → NPT → Production pipeline with auto-generated run scripts and GPU support |
| **Analysis** | RMSD, RMSF, contacts, H-bonds, SASA, clustering, PCA, interactive HTML reports |
| **PMF** | Steered MD setup, automated umbrella sampling window generation, WHAM analysis |
| **MM/PBSA** | Single-frame and trajectory binding energy with electrostatic, vdW, solvation decomposition |

## System Requirements

!!! info "Required"
    - **Python** 3.8 - 3.11 (3.10 recommended)
    - **GROMACS** 2024.3+ (CUDA support recommended)
    - **PDBFixer**, **NumPy**, **SciPy**

!!! tip "Force Field Specific"
    - **GAFF/GAFF2**: AmberTools, ACPYPE
    - **OpenFF**: OpenFF Toolkit, OpenFF Interchange
    - **CGenFF**: Downloaded parameters from ParamChem
    - **OPLS-AA**: Internet connection (LigParGen)

## Documentation

- **[Getting Started](getting-started/index.md)** -- Installation and basic usage
- **[User Guide](user-guide/index.md)** -- Comprehensive feature documentation
- **[Tutorials](tutorials/index.md)** -- Step-by-step guides and examples
- **[API Reference](api/index.md)** -- Complete API documentation
- **[Examples](examples/index.md)** -- Real-world use cases and scripts

## Contributing

```bash
git clone https://github.com/your-username/PRISM.git
cd PRISM && pip install -e .[dev]
pytest tests/
```

For questions or suggestions, open an issue on [GitHub](https://github.com/AIB001/PRISM/issues).

## Citation

```bibtex
@software{prism2024,
  author = {Shi, Zhaoqi},
  title = {PRISM: Protein Receptor Interaction Simulation Modeler},
  year = {2024},
  version = {1.2.0},
  institution = {Institute of Quantitative Biology, Zhejiang University},
  url = {https://github.com/AIB001/PRISM}
}
```

---

<div style="text-align: center; opacity: 0.7;" markdown>

**Institute of Quantitative Biology** | Zhejiang University | **Author:** Zhaoqi Shi

</div>
