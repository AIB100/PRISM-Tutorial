# Changelog

All notable changes to PRISM will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [1.2.0] - 2024-10-25

### Added
- **PMF Module**: Complete umbrella sampling and WHAM analysis workflow
  - Automated SMD generation
  - Umbrella window setup and execution
  - WHAM analysis with bootstrap error estimation
  - Binding free energy calculations
- **Batch Processing**: High-throughput virtual screening support
  - Parallel system building
  - Batch trajectory processing
  - Multi-system analysis and comparison
- **Interactive Visualization**: HTML-based contact analysis
  - 2D/3D interactive contact maps
  - Drag-and-drop interface
  - High-resolution image export
- **Trajectory Processing**: PBC artifact correction
  - DCD to XTC conversion
  - Automatic centering and wrapping
  - Batch processing support
- **Force Field Support**: Added MMFF, MATCH, and Hybrid options
- **Configuration System**: YAML-based configuration files
  - Template export functionality
  - Parameter validation
  - Configuration inheritance

### Changed
- **Improved Protein Preparation**: Enhanced metal ion handling
  - Smart mode for selective metal retention
  - Distance-based filtering
  - Better artifact removal
- **Enhanced Error Messages**: More informative error reporting
- **Updated Documentation**: Comprehensive tutorials and examples
- **Performance**: Optimized ligand parameterization speed

### Fixed
- Issue with OPLS-AA force field generation
- Box size calculation for large systems
- Ion concentration accuracy
- Trajectory analysis memory usage

---

## [1.1.0] - 2024-08-15

### Added
- **Analysis Module**: Comprehensive trajectory analysis tools
  - RMSD/RMSF calculations
  - Contact frequency analysis
  - Hydrogen bond detection
  - SASA calculations
- **Multiple Force Fields**: Support for additional ligand force fields
  - OpenFF (Open Force Field)
  - OPLS-AA (via LigParGen)
  - CGenFF (CHARMM)
- **Command-Line Interface**: Full CLI for all operations
- **Progress Indicators**: Real-time build progress display

### Changed
- Refactored builder architecture for modularity
- Improved GROMACS version detection
- Enhanced water model selection

### Fixed
- GAFF2 parameter generation edge cases
- Solvation box shape issues
- MDP file generation for CHARMM force fields

---

## [1.0.0] - 2024-06-01

### Added
- **Initial Release**: Core system building functionality
  - Protein-ligand system construction
  - GAFF/GAFF2 force field support
  - AMBER protein force fields
  - Automated solvation and neutralization
  - MDP file generation
  - GROMACS integration
- **Documentation**: User guide and API reference
- **Examples**: Basic tutorials and examples

---

## [0.9.0] - 2024-04-15 (Beta)

### Added
- Beta testing version
- Core PRISMBuilder class
- Basic force field generation
- Protein cleaning utilities

### Known Issues
- Limited force field options
- No analysis tools
- Documentation incomplete

---

## Planned Features

### [1.3.0] - Planned Q1 2025
- **MM/PBSA Module**: Binding energy calculations
- **Mutation Analysis**: Automated mutation scanning
- **Enhanced Visualization**: Advanced 3D visualization
- **Machine Learning Integration**: ML-based predictions
- **Cloud Computing**: AWS/GCP integration

### [2.0.0] - Planned Q3 2025
- **Multi-Ligand Support**: Handle multiple ligands
- **Membrane Proteins**: Automated membrane system building
- **Alchemical Free Energy**: FEP/TI support
- **Enhanced PMF**: Multiple reaction coordinates
- **Web Interface**: Browser-based interface

---

## Version History Summary

| Version | Release Date | Major Features |
|---------|--------------|----------------|
| 1.2.0 | 2024-10-25 | PMF, Batch Processing, Visualization |
| 1.1.0 | 2024-08-15 | Analysis Tools, More Force Fields |
| 1.0.0 | 2024-06-01 | Initial Release |
| 0.9.0 | 2024-04-15 | Beta Version |

---

## Upgrade Guide

### From 1.1.0 to 1.2.0

**Breaking Changes**: None

**New Features**:
```python
# New PMF functionality
import prism
results = prism.run_pmf_workflow("md_dir", "pmf_output")

# New batch processing
prism.batch_process_trajectories(trajectories, "output", "topology.tpr")

# New visualization
prism.visualize_trajectory("traj.xtc", "system.gro", "ligand.sdf")
```

**Deprecated**: Nothing deprecated

### From 1.0.0 to 1.1.0

**Breaking Changes**: None

**New Features**:
```python
# Analysis module
analyzer = prism.analyze_trajectory("system.gro", "traj.xtc")

# Additional force fields
system = prism.system(
    "protein.pdb", "ligand.sdf",
    ligand_forcefield="openff"  # New option
)
```

---

## Migration Notes

### Configuration Files

1.2.0 introduces YAML configuration files. Old command-line arguments still work:

```python
# Old style (still works)
system = prism.system("protein.pdb", "ligand.mol2",
                      forcefield="amber99sb",
                      water_model="tip3p")

# New style (recommended)
system = prism.system("protein.pdb", "ligand.mol2",
                      config="my_config.yaml")
```

### Analysis

The analysis module API is stable but new functions are added:

```python
# 1.1.0
analyzer.calc_rmsd()
analyzer.calc_contacts()

# 1.2.0 (added)
analyzer.calc_dihedral()
analyzer.calc_clustering()
```

---

## Reporting Issues

Found a bug or regression?

1. Check if it's already reported: [GitHub Issues](https://github.com/AIB001/PRISM/issues)
2. Search discussions: [GitHub Discussions](https://github.com/AIB001/PRISM/discussions)
3. Create new issue with:
   - PRISM version (`prism.get_version()`)
   - Python version
   - Operating system
   - Minimal example to reproduce

---

## Contributing

For contributions, please open an issue or pull request on [GitHub](https://github.com/AIB001/PRISM).

---

**Full Changelog**: [github.com/AIB001/PRISM/blob/main/CHANGELOG.md](https://github.com/AIB001/PRISM/blob/main/CHANGELOG.md)

**Last Updated**: October 25, 2024
