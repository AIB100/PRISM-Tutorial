# License

PRISM is released under the MIT License.

## MIT License

```
MIT License

Copyright (c) 2024 Zhaoqi Shi
Theoretical Chemistry Institute
University of Wisconsin-Madison

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## What This Means

### You CAN:

- ✓ Use PRISM for **any purpose** (commercial or academic)
- ✓ **Modify** the source code
- ✓ **Distribute** original or modified versions
- ✓ **Incorporate** into proprietary software
- ✓ Use in **closed-source** projects

### You MUST:

- Include the **copyright notice** in copies
- Include the **license text** in distributions

### You CANNOT:

- Hold the authors **liable** for damages
- Claim the software comes with a **warranty**

---

## Dependencies

PRISM uses several third-party tools and libraries, each with their own licenses:

### Required Dependencies

**GROMACS**
- License: LGPL 2.1
- Website: [gromacs.org](https://www.gromacs.org)
- Use: MD simulation engine

**Python Libraries**
- NumPy: BSD License
- SciPy: BSD License
- Matplotlib: BSD-compatible

### Optional Dependencies

**AmberTools (for GAFF/GAFF2)**
- License: GPL/LGPL
- Website: [ambermd.org](https://ambermd.org)

**OpenFF Toolkit (for OpenFF)**
- License: MIT
- Website: [openforcefield.org](https://openforcefield.org)

**MDTraj (for analysis)**
- License: LGPL 2.1
- Website: [mdtraj.org](https://www.mdtraj.org)

**RDKit (for visualization)**
- License: BSD
- Website: [rdkit.org](https://www.rdkit.org)

---

## Force Field Licenses

### GAFF / GAFF2
- Distributed with AmberTools
- Free for academic and commercial use
- License: GPL/LGPL

### OpenFF
- Open Force Field Initiative
- License: MIT
- Free for all uses

### OPLS-AA (via LigParGen)
- Web service provided by Jorgensen Lab
- Free for academic use
- Commercial use: Contact developers

### CGenFF
- CHARMM General Force Field
- Free download for academic use
- Commercial use: License required

---

## Citation Requirement

While not legally required by the MIT License, we kindly request that you cite PRISM in any publications or presentations that use the software:

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

See [Publications](publications.md) for details on citing PRISM.

---

## Contributing

By contributing to PRISM, you agree that your contributions will be licensed under the same MIT License.

See [Contributing Guide](../development/contributing.md) for details.

---

## Disclaimer

PRISM is research software provided "as is" without warranty. While we strive for accuracy and reliability:

- No guarantee of correctness for scientific results
- No guarantee of fitness for any particular purpose
- Use at your own risk
- Validate results independently

Always:
- Test on known systems
- Compare with experimental data when available
- Use multiple methods for critical results
- Understand the underlying science

---

## Contact

Questions about licensing?

- **Email**: zshi268@wisc.edu
- **GitHub**: [github.com/AIB001/PRISM](https://github.com/AIB001/PRISM)

---

## Other Licenses

PRISM includes or may use code from:

### ACPYPE
- License: GPL v3
- Used for: GAFF parameter conversion
- Website: [github.com/alanwilter/acpype](https://github.com/alanwilter/acpype)

### PDBFixer
- License: MIT
- Used for: Protein structure preparation
- Website: [github.com/openmm/pdbfixer](https://github.com/openmm/pdbfixer)

All such code is used in compliance with their respective licenses.

---

**Full License**: [github.com/AIB001/PRISM/blob/main/LICENSE](https://github.com/AIB001/PRISM/blob/main/LICENSE)

**Last Updated**: October 2024
