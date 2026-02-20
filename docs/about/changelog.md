# Changelog

All notable changes to PRISM will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [2.0.0] - 2025-02

### Added
- **PMF Module**: Complete umbrella sampling and WHAM analysis workflow
  - Metropolis-Hastings optimized pulling direction selection
  - Automated SMD, umbrella window setup, and WHAM analysis
  - Binding free energy calculations with bootstrap error estimation
- **REST2 Module**: Replica Exchange with Solute Tempering 2
  - Geometric temperature ladder generation
  - Automatic hot region identification and topology scaling
  - Multi-replica Hamiltonian exchange setup
- **MM/PBSA Module**: Binding energy calculations
  - Single-frame and trajectory-based modes
  - gmx_MMPBSA and AMBER MMPBSA.py backends
  - Automatic GROMACS-to-AMBER topology conversion
- **Gaussian RESP Charges**: HF/6-31G* and DFT charge calculations
  - Gaussian input file generation
  - Two-stage RESP fitting
  - Integration with GAFF/GAFF2 parameterization
- **Batch Processing**: High-throughput virtual screening support
- **Interactive Visualization**: HTML-based contact analysis
- **Trajectory Processing**: PBC artifact correction and format conversion
- **Multi-Ligand Support**: Handle multiple ligands in a single system
- **Configuration System**: YAML-based configuration files

### Changed
- Improved protein preparation with enhanced metal ion handling
- Enhanced error messages with colored terminal output
- Optimized ligand parameterization speed
- Improved meeko dependency handling and protonation mode

### Fixed
- OPLS-AA force field generation issues
- Box size calculation for large systems
- Ion concentration accuracy
- Trajectory analysis memory usage

---

## Reporting Issues

Found a bug?

1. Check if it's already reported: [GitHub Issues](https://github.com/AIB001/PRISM/issues)
2. Create new issue with:
   - PRISM version (`prism.get_version()`)
   - Python version
   - Operating system
   - Minimal example to reproduce

---

## Contributing

For contributions, please open an issue or pull request on [GitHub](https://github.com/AIB001/PRISM).

---

**Last Updated**: February 2025
