# Examples Overview

This section provides practical, real-world examples demonstrating various PRISM capabilities. Each example includes complete, runnable code and detailed explanations.

## Available Examples

### Basic Examples

<div class="grid cards" markdown>

-   **Simple Protein-Ligand**

    ---

    Build and simulate a basic protein-ligand complex

    [Go to example](simple.md)

-   **PMF Calculation**

    ---

    Estimate binding free energy via steered MD and umbrella sampling

    [Go to example](pmf-example.md)

-   **Multiple Force Fields**

    ---

    Compare results across different force fields

    [Go to example](multi-forcefield.md)

-   **Custom Configuration**

    ---

    Advanced system setup with custom parameters

    [Go to example](custom-config.md)

</div>

---

## Example Structure

Each example follows this structure:

1. **Overview** - What the example demonstrates
2. **Prerequisites** - Required files and dependencies
3. **Complete Code** - Full working script
4. **Step-by-Step Explanation** - Detailed walkthrough
5. **Expected Results** - What you should see
6. **Troubleshooting** - Common issues and solutions

---

## Quick Start Template

Use this template as a starting point for your own projects:

```python
#!/usr/bin/env python3
"""
PRISM Example Template
"""

import prism

# 1. Define inputs
protein_file = "protein.pdb"
ligand_file = "ligand.mol2"
output_directory = "my_simulation"

# 2. Create system
system = prism.system(
    protein_file,
    ligand_file,
    output_dir=output_directory,
    ligand_forcefield="gaff",
    forcefield="amber99sb",
    water_model="tip3p"
)

# 3. Build system
print("Building system...")
system.build()

# 4. Display information
print("\nSystem built successfully!")
print(f"Output directory: {output_directory}")
print(f"Next step: cd {output_directory}/GMX_PROLIG_MD && bash localrun.sh")
```

---

## Example Datasets

All examples use freely available test systems:

### T4 Lysozyme - Benzene
- **Protein**: T4 Lysozyme L99A/M102Q (PDB: [5JWT](https://www.rcsb.org/structure/5JWT))
- **Ligand**: Benzene
- **Use case**: Small molecule binding, hydrophobic interactions
- **Simulation time**: ~2-4 hours on GPU

### CDK2 - CVT-313
- **Protein**: Cyclin-dependent kinase 2 (PDB: [6INL](https://www.rcsb.org/structure/6INL), Swiss-Model homology model)
- **Ligand**: CVT-313 (purine-based inhibitor)
- **Use case**: PMF binding free energy, Gaussian RESP charges
- **Simulation time**: days (umbrella sampling, parallelizable)

### BACE1 - Inhibitor
- **Protein**: Beta-secretase 1 (PDB: 3L5D)
- **Ligand**: Small molecule inhibitor
- **Use case**: Drug design, hydrogen bonding
- **Simulation time**: ~4-6 hours on GPU

### HSP90 - Geldanamycin
- **Protein**: Heat Shock Protein 90 (PDB: 1YET)
- **Ligand**: Geldanamycin
- **Use case**: Large ligand, complex binding
- **Simulation time**: ~6-8 hours on GPU

---

## Running Examples

### Option 1: Download and Run

```bash
# Clone examples repository
git clone https://github.com/AIB001/PRISM-examples.git
cd PRISM-examples

# Run specific example
cd simple-protein-ligand
python run_example.py
```

### Option 2: Copy-Paste Code

Each example page includes complete, copy-paste ready code that you can run directly.

### Option 3: Jupyter Notebooks

Interactive Jupyter notebooks are available for all examples:

```bash
cd PRISM-examples
jupyter notebook
```

---

## Example Categories

### System Building
- [Simple Protein-Ligand](simple.md) - Basic workflow
- [Custom Configuration](custom-config.md) - Advanced setup
- [Multiple Force Fields](multi-forcefield.md) - Force field comparison

### Analysis
- Trajectory analysis and visualization
- Contact frequency mapping
- RMSD/RMSF calculations
- Hydrogen bond analysis

### Advanced Calculations
- [PMF Calculation](pmf-example.md) - Binding free energy via umbrella sampling and WHAM
- MM/PBSA binding energies
- Free energy perturbation

### High-Throughput
- Batch processing
- Virtual screening workflows
- Parallel execution

---

## Tips for Using Examples

1. **Start Simple**: Begin with the [Simple Example](simple.md)
2. **Modify Incrementally**: Change one parameter at a time
3. **Check Dependencies**: Ensure all required tools are installed
4. **Use Test Systems**: Start with small systems before scaling up
5. **Read Error Messages**: They usually point to the issue

---

## Contributing Examples

Have a useful example to share? We welcome contributions!

1. Fork the PRISM repository
2. Create your example in `examples/`
3. Follow the standard structure
4. Submit a pull request

For contributions, please open an issue or pull request on [GitHub](https://github.com/AIB001/PRISM).

---

## Getting Help

If you encounter issues with examples:

1. Check the [Troubleshooting Guide](../user-guide/troubleshooting.md)
2. Search [GitHub Issues](https://github.com/AIB001/PRISM/issues)
3. Ask in [GitHub Discussions](https://github.com/AIB001/PRISM/discussions)
4. Email: zshi268@wisc.edu

---

<div class="whats-next" markdown>

## What's Next

- [Try the Simple Example to get started](simple.md)
- [Explore Tutorials for guided learning](../tutorials/index.md)
- [Check the API Reference for scripting details](../api/index.md)

</div>
