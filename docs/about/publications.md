# Publications

Research using PRISM and related publications.

## PRISM Software

### Primary Citation

If you use PRISM in your research, please cite:

```bibtex
@software{prism2024,
  author = {Shi, Zhaoqi},
  title = {PRISM: Protein Receptor Interaction Simulation Modeler},
  year = {2024},
  version = {2.0.0},
  institution = {Institute of Quantitative Biology, Zhejiang University},
  url = {https://github.com/AIB001/PRISM}
}
```

### Related Publications

Publications describing PRISM methodology and applications (in preparation):

1. **PRISM: An Automated Toolkit for Protein-Ligand Molecular Dynamics Simulations**  
   Shi, Z. et al.  
   *Manuscript in preparation*, 2024

2. **Force Field Comparison for Protein-Ligand Binding: Insights from PRISM**  
   Shi, Z. et al.  
   *Manuscript in preparation*, 2024

---

## Methodology Papers

### Force Fields

**GAFF (General AMBER Force Field)**
- Wang, J., et al. (2004). Development and testing of a general amber force field. *J. Comput. Chem.*, 25(9), 1157-1174.
- DOI: [10.1002/jcc.20035](https://doi.org/10.1002/jcc.20035)

**GAFF2**
- He, X., et al. (2020). A fast and high-quality charge model for the next generation general AMBER force field. *J. Chem. Phys.*, 153(11), 114502.

**OpenFF (Open Force Field)**
- Mobley, D. L., et al. (2020). Open Force Field Initiative. *Living J. Comp. Mol. Sci.*, 2(1), 18378.

**OPLS-AA**
- Jorgensen, W. L., et al. (1996). Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids. *J. Am. Chem. Soc.*, 118(45), 11225-11236.

### Free Energy Methods

**Umbrella Sampling**
- Torrie, G. M., & Valleau, J. P. (1977). Nonphysical sampling distributions in Monte Carlo free-energy estimation: Umbrella sampling. *J. Comput. Phys.*, 23(2), 187-199.

**WHAM**
- Kumar, S., et al. (1992). The weighted histogram analysis method for free-energy calculations on biomolecules. I. The method. *J. Comput. Chem.*, 13(8), 1011-1021.
- Hub, J. S., et al. (2010). g_wham—A Free Weighted Histogram Analysis Implementation Including Robust Error and Autocorrelation Estimates. *J. Chem. Theory Comput.*, 6(12), 3713-3720.

**PMF Calculations**
- Lemkul, J. A. (2019). From Proteins to Perturbed Hamiltonians: A Suite of Tutorials for the GROMACS-2018 Molecular Simulation Package. *Living J. Comp. Mol. Sci.*, 1(1), 5068.

### MD Simulations

**GROMACS**
- Abraham, M. J., et al. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. *SoftwareX*, 1-2, 19-25.
- Van Der Spoel, D., et al. (2005). GROMACS: fast, flexible, and free. *J. Comput. Chem.*, 26(16), 1701-1718.

---

## Submit Your Publication

If you use PRISM in your research, we'd be happy to list your publication here.

- Email: zshi268@wisc.edu
- Subject: "PRISM Publication"
- Include: Title, Authors, Journal, DOI

---

## Citation Guidelines

### In Methods Section

Example text for methods section:

> "Protein-ligand systems were built using PRISM (Protein Receptor Interaction Simulation Modeler) version 2.0.0 [citation]. Ligand force field parameters were generated using [GAFF/OpenFF/etc], and the protein was modeled with the AMBER99SB force field. The system was solvated in a TIP3P water box with 0.15 M NaCl. Molecular dynamics simulations were performed using GROMACS 2024.3."

### In Acknowledgments

> "We thank Zhaoqi Shi for developing PRISM and providing support."

---

## Related Software

### Dependencies
- **GROMACS**: [gromacs.org](https://www.gromacs.org)
- **AmberTools**: [ambermd.org](https://ambermd.org)
- **OpenFF**: [openforcefield.org](https://openforcefield.org)
- **MDTraj**: [mdtraj.org](https://www.mdtraj.org)

### Complementary Tools
- **VMD**: Visualization
- **PyMOL**: Molecular graphics
- **ProDy**: Protein dynamics analysis
- **MDAnalysis**: Trajectory analysis

---

**Last Updated**: February 2025
