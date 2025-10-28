# Welcome to PRISM Documentation

<p align="center">
  <img src="assets/logo.png" alt="PRISM Logo" width="300">
</p>

<p align="center">
  <strong>Protein Receptor Interaction Simulation Modeler</strong><br>
  <em>A comprehensive toolkit for building and analyzing protein-ligand MD systems</em>
</p>

<p align="center">
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/Python-3.8%2B-blue.svg" alt="Python Version"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License"></a>
  <a href="https://prism-docs.github.io"><img src="https://img.shields.io/badge/Docs-Latest-brightgreen.svg" alt="Documentation"></a>
  <a href="https://github.com/AIB001/PRISM"><img src="https://img.shields.io/badge/GitHub-PRISM-black.svg" alt="GitHub"></a>
</p>

---

## Overview

**PRISM** (Protein Receptor Interaction Simulation Modeler) is a powerful Python framework designed for **building protein-ligand systems** and **analyzing MD trajectories**, developed at Theoretical Chemistry Institute, University of Wisconsin-Madison. It automates the complex workflow of preparing molecular dynamics simulations with multiple force field options and provides comprehensive analysis tools for trajectory data.

## Quick Start

### Installation

Install from source (recommended):

```bash
git clone https://github.com/AIB001/PRISM.git
cd PRISM
pip install -e .
```

### Basic Usage

Build a protein-ligand system with just one command:

```python
import prism

# Build system with GAFF force field (default)
system = prism.system("protein.pdb", "ligand.mol2", output_dir="my_system")
system.build()
```

Or use the command-line interface:

```bash
# Using GAFF force field (default)
prism protein.pdb ligand.mol2 -o my_system

# Using OpenFF force field
prism protein.pdb ligand.sdf -o my_system --ligand-forcefield openff

# Using GAFF2 (improved GAFF)
prism protein.pdb ligand.mol2 -o my_system --ligand-forcefield gaff2
```

### Run MD Simulation

After building, run the simulation:

```bash
cd my_system/GMX_PROLIG_MD
bash localrun.sh
```

## Key Features

<div class="feature-grid">
  <div class="feature-card">
    <h3>🔧 Automated System Building</h3>
    <p>One-command setup for protein-ligand MD systems with automated protein cleaning, ligand parameterization, solvation, and neutralization.</p>
  </div>

  <div class="feature-card">
    <h3>🧪 Multiple Force Fields</h3>
    <p>Support for 8+ ligand force fields: GAFF, GAFF2, OpenFF, CGenFF, OPLS-AA, MMFF, MATCH, and hybrid approaches.</p>
  </div>

  <div class="feature-card">
    <h3>📈 Trajectory Analysis</h3>
    <p>Comprehensive analysis tools including RMSD/RMSF, contacts, hydrogen bonds, clustering, SASA, and dihedral angles.</p>
  </div>

  <div class="feature-card">
    <h3>📊 Interactive Visualization</h3>
    <p>Generate interactive HTML visualizations of protein-ligand contacts with customizable 2D/3D views.</p>
  </div>

  <div class="feature-card">
    <h3>⚡ GROMACS Integration</h3>
    <p>Native GROMACS workflow with auto-generated MDP files for EM, NVT, NPT, and production MD.</p>
  </div>

  <div class="feature-card">
    <h3>🎯 Advanced Calculations</h3>
    <p>Built-in support for PMF (umbrella sampling/WHAM) and MM/PBSA binding energy calculations.</p>
  </div>
</div>

## Core Modules

### Builder Module
Automated system construction for protein-ligand complexes:

- **Force Field Generation**: Automatic ligand parameterization with GAFF, GAFF2, OpenFF, CGenFF, OPLS-AA, etc.
- **Protein Preparation**: Intelligent cleaning with smart metal ion handling, protonation optimization
- **System Assembly**: Solvation, neutralization, ion addition with configurable box shapes
- **GROMACS Integration**: Direct topology and coordinate file generation

### Simulation Module
Interface for running and managing MD simulations:

- **Workflow Automation**: Complete EM → NVT → NPT → Production pipeline
- **Auto-generated Scripts**: Ready-to-run bash scripts with checkpoint restart support
- **GROMACS Integration**: Native support with GPU acceleration options

### Analysis Module
Comprehensive trajectory analysis tools:

- **Geometric Analysis**: RMSD, RMSF, distance calculations
- **Structural Analysis**: Contacts, hydrogen bonds, SASA, dihedral angles
- **Advanced Analysis**: Clustering, PCA, correlation analysis
- **Visualization**: Interactive HTML reports with 2D/3D contact maps

### PMF Module
Binding free energy calculations:

- **SMD Preparation**: Steered molecular dynamics setup
- **Umbrella Sampling**: Automated window generation and equilibration
- **WHAM Analysis**: Potential of mean force calculation
- **Energy Profiles**: Binding energy estimation

### MM/PBSA Module
Molecular mechanics Poisson-Boltzmann surface area:

- **Single-frame**: Docking pose energy evaluation
- **Trajectory**: Time-averaged binding energy from MD
- **Component Analysis**: Decomposition into electrostatic, vdW, polar/non-polar solvation

## System Requirements

### Required Dependencies
- **Python**: 3.8 - 3.11 (3.10 recommended)
- **GROMACS**: 2024.3 or later (with CUDA support recommended)
- **PDBFixer**: For protein structure preparation
- **NumPy, SciPy**: For numerical computations

### Force Field Specific
- **GAFF/GAFF2**: AmberTools, ACPYPE
- **OpenFF**: OpenFF Toolkit, OpenFF Interchange
- **CGenFF**: Downloaded CGenFF parameters from ParamChem
- **OPLS-AA**: Internet connection (uses LigParGen web server)
- **MMFF/MATCH**: Internet connection (uses SwissParam web server)

### Optional for Analysis
- **MDTraj**: Trajectory analysis and visualization
- **RDKit**: Enhanced molecular structure handling
- **Matplotlib/Seaborn**: Data visualization

### Hardware
- **Memory**: 8 GB RAM minimum (16 GB+ recommended for large systems)
- **Storage**: 10+ GB for typical protein-ligand system
- **GPU**: NVIDIA GPU with CUDA support (highly recommended for MD simulations)

## Documentation Structure

- **[Getting Started](getting-started/index.md)**: Installation and basic usage
- **[User Guide](user-guide/index.md)**: Comprehensive feature documentation
- **[Tutorials](tutorials/index.md)**: Step-by-step guides and examples
- **[API Reference](api/index.md)**: Complete API documentation
- **[Examples](examples/index.md)**: Real-world use cases and scripts

## Contributing

We welcome contributions! To get started:

```bash
# Fork and clone repository
git clone https://github.com/your-username/PRISM.git
cd PRISM

# Install in development mode
pip install -e .[dev]

# Run tests
pytest tests/

# Submit pull request via GitHub
```

For questions or suggestions, please open an issue on [GitHub](https://github.com/AIB001/PRISM/issues).

## Citation

If you use PRISM in your research, please cite:

```bibtex
@software{prism2024,
  author = {Shi, Zhaoqi},
  title = {PRISM: Protein Receptor Interaction Simulation Modeler},
  year = {2024},
  version = {1.2.0},
  institution = {Theoretical Chemistry Institute, University of Wisconsin-Madison},
  url = {https://github.com/AIB001/PRISM}
}
```

## Support

- **Documentation**: [https://prism-docs.github.io](https://prism-docs.github.io)
- **Issue Tracker**: [GitHub Issues](https://github.com/AIB001/PRISM/issues)
- **Discussions**: [GitHub Discussions](https://github.com/AIB001/PRISM/discussions)
- **Email**: zshi268@wisc.edu

## License

PRISM is released under the [MIT License](about/license.md).

---

<p align="center">
  <strong>Developed at</strong><br>
  Theoretical Chemistry Institute<br>
  University of Wisconsin-Madison<br>
  <br>
  <strong>Author:</strong> Zhaoqi Shi
</p>