# About PRISM

**PRISM** (Protein Receptor Interaction Simulation Modeler) is a comprehensive Python toolkit for building and analyzing protein-ligand molecular dynamics systems.

## Mission

PRISM's mission is to make high-quality molecular dynamics simulations accessible to researchers by automating the complex workflow of system preparation, simulation setup, and trajectory analysis.

## History

PRISM was developed at the Institute of Quantitative Biology, Zhejiang University, starting in 2023 to address the need for:

- **Automation**: Eliminate manual, error-prone system building steps
- **Flexibility**: Support multiple force field options
- **Reproducibility**: Consistent, documented workflows
- **Accessibility**: Lower barriers for non-experts

## Key Features

### System Building
- Automated protein-ligand system construction
- 8+ ligand force field options (GAFF, GAFF2, OpenFF, OPLS-AA, CGenFF, MMFF, MATCH, Hybrid)
- Multiple protein force fields (AMBER, CHARMM, OPLS)
- Intelligent protein preparation with metal ion handling
- One-command GROMACS integration

### Analysis
- Comprehensive trajectory analysis (RMSD, RMSF, contacts, H-bonds)
- Interactive HTML visualizations
- Automated plotting and reporting
- Publication-quality figures

### Advanced Calculations
- PMF/Umbrella sampling for binding free energies
- MM/PBSA calculations
- Batch processing for high-throughput screening
- Parallel execution support

## Design Philosophy

1. **Simple for Simple Tasks**: Common workflows should be one command
2. **Flexible for Complex Tasks**: Full control when needed
3. **Reproducible**: All settings documented and version controlled
4. **Well-Documented**: Clear examples and comprehensive docs
5. **Open Source**: Free and community-driven

## Technology Stack

- **Language**: Python 3.8-3.11
- **MD Engine**: GROMACS 2024.3+
- **Force Fields**: AmberTools, OpenFF, LigParGen, ParamChem
- **Analysis**: MDTraj, NumPy, SciPy, Matplotlib
- **Visualization**: RDKit, Plotly, custom HTML generators

## Use Cases

### Drug Discovery
- Virtual screening of compound libraries
- Binding affinity prediction
- Force field validation

### Structural Biology
- Protein-ligand binding mechanisms
- Conformational changes
- Allosteric effects

### Method Development
- Force field comparison and benchmarking
- Protocol optimization
- Algorithm testing

## Community

PRISM is developed and maintained by researchers, for researchers. We welcome:

- **Bug Reports**: Help us improve reliability
- **Feature Requests**: Tell us what you need
- **Contributions**: Code, docs, examples
- **Feedback**: What works, what doesn't

## Getting Involved

- **GitHub**: [github.com/AIB001/PRISM](https://github.com/AIB001/PRISM)
- **Discussions**: Share ideas and get help
- **Issues**: Report bugs and request features
- **Email**: Contact the development team

## Citing PRISM

If you use PRISM in your research, please cite:

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

## Acknowledgments

PRISM development has been supported by:
- Institute of Quantitative Biology, Zhejiang University
- Open source community contributions
- Feedback from users worldwide

## License

PRISM is released under the MIT License, allowing free use in academic and commercial settings. See [License](license.md) for details.

---

**Developed at:**
Institute of Quantitative Biology
Zhejiang University

**Principal Developer:**  
Zhaoqi Shi

**Contact:**  
zshi268@wisc.edu

