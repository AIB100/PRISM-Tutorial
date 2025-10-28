# Installation Guide

This guide provides detailed instructions for installing PRISM and its dependencies on various platforms.

## 📋 Prerequisites

PRISM requires several components to be installed. Install them in the following order:

### 1. GROMACS (Required)

GROMACS is essential for PRISM to build and run molecular dynamics simulations. PRISM uses GROMACS for:
- Processing protein structures (`pdb2gmx`)
- Building simulation boxes and solvation (`editconf`, `solvate`)
- Adding ions (`genion`)
- Running MD simulations (`mdrun`)

Follow these steps to install GROMACS with GPU support:

#### Check CUDA Toolkit
First, verify that you have CUDA installed:

```bash
nvcc --version
```

If CUDA is not installed, please install it from [NVIDIA's CUDA Toolkit page](https://developer.nvidia.com/cuda-downloads) before proceeding.

#### Install Build Dependencies

=== "Ubuntu/Debian"

    ```bash
    # Install essential build tools
    sudo apt-get update
    sudo apt-get install gcc g++ cmake
    sudo apt-get install build-essential
    ```

=== "CentOS/RHEL"

    ```bash
    # Install essential build tools
    sudo yum groupinstall "Development Tools"
    sudo yum install cmake gcc-c++
    ```

=== "macOS"

    ```bash
    # Install using Homebrew
    brew install cmake gcc
    ```

#### Download and Build GROMACS

```bash
# Download GROMACS 2024.3
wget https://ftp.gromacs.org/gromacs/gromacs-2024.3.tar.gz
tar xfz gromacs-2024.3.tar.gz

# Create build directory
cd gromacs-2024.3
mkdir build
cd build

# Configure with CMake (adjust paths as needed)
cmake .. -DGMX_MPI=ON \
         -DGMX_BUILD_OWN_FFTW=ON \
         -DGMX_GPU=CUDA \
         -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
         -DCUDA_INCLUDE_DIRS=/usr/local/cuda/include \
         -DCUDA_CUDART_LIBRARY=/usr/local/cuda/lib64 \
         -DCMAKE_INSTALL_PREFIX=/opt/gromacs-2024.3

# Build GROMACS (use all available cores)
make -j$(nproc)

# Run tests (optional but recommended)
make check

# Install GROMACS
sudo make install

# Source GROMACS environment
source /opt/gromacs-2024.3/bin/GMXRC
```

!!! tip "Add to your shell configuration"
    Add the following line to your `~/.bashrc` or `~/.zshrc` to automatically source GROMACS:
    ```bash
    source /opt/gromacs-2024.3/bin/GMXRC
    ```

#### Install Bioconda (Optional)

For additional bioinformatics tools:

```bash
# Install bioconda channel
conda config --add channels bioconda
conda config --add channels conda-forge
```

### 2. Python Environment (Required)

PRISM requires Python 3.10 or higher. We recommend using Conda for environment management:

```bash
# Create a new environment for PRISM
conda create -n prism python=3.10
conda activate prism
```

!!! info "Python Version"
    While Python 3.10 is recommended, PRISM supports Python 3.8 through 3.11. Python 3.10 offers the best balance of compatibility and performance.

### 3. Core Python Dependencies (Required)

Install essential Python packages:

```bash
# Activate PRISM environment
conda activate prism

# Install PDBFixer and core scientific packages
conda install -c conda-forge pdbfixer numpy scipy

# Install YAML support
pip install pyyaml
```

## 🧪 Force Field Specific Dependencies

PRISM supports **8+ ligand force fields**. Install dependencies based on your needs:

### For GAFF / GAFF2 Support (Most Common)

The General AMBER Force Field is the default and most widely used:

```bash
# AmberTools (provides antechamber, parmchk2)
conda install -c conda-forge ambertools

# ACPYPE for GROMACS conversion (required)
pip install acpype

# RDKit for molecular structure handling (recommended)
conda install -c conda-forge rdkit
```

### For OpenFF Support

Modern, systematically improved force fields:

```bash
# OpenFF toolkit and interchange
conda install -c conda-forge openff-toolkit openff-interchange

# RDKit for molecular reading
conda install -c conda-forge rdkit
```

### For CGenFF Support

CHARMM General Force Field (requires web download):

1. Upload your ligand to [ParamChem](https://cgenff.umaryland.edu/)
2. Download the generated parameter files
3. Provide the directory path to PRISM with `--forcefield-path`

```bash
# No additional software needed - uses downloaded files
prism protein.pdb ligand.mol2 -o output --ligand-forcefield cgenff --forcefield-path /path/to/cgenff_files
```

### For OPLS-AA Support

Uses LigParGen web server (requires internet):

```bash
# RDKit for structure handling
conda install -c conda-forge rdkit

# No additional software - uses LigParGen API
```

### For MMFF / MATCH / Hybrid Support

Uses SwissParam web server (requires internet):

```bash
# RDKit for structure handling
conda install -c conda-forge rdkit

# No additional software - uses SwissParam API
```

## 📦 Installing PRISM

Once all prerequisites are installed, you can install PRISM:

### Recommended: Development Installation

This allows you to modify PRISM and immediately see changes:

```bash
# Clone the PRISM repository
git clone https://github.com/AIB001/PRISM.git
cd PRISM

# Install in development mode
pip install -e .

# Verify installation
prism --help
python -c "import prism; print(prism.__version__)"
```

### Alternative: Direct Usage (Without Installation)

You can also use PRISM directly:

```bash
# Clone the repository
git clone https://github.com/AIB001/PRISM.git
cd PRISM

# Run PRISM module directly
python -m prism.builder protein.pdb ligand.mol2 -o output_dir
```

## ✅ Verification

After installation, verify that everything is working correctly:

### Test GROMACS Installation

```bash
# Check GROMACS version
gmx --version

# Should show version 2024.3 with CUDA support
```

### Test Python Environment

```bash
# Activate PRISM environment
conda activate prism

# Test core functionality
python -c "import prism; print(f'PRISM version: {prism.__version__}')"
python -c "import pdbfixer; print('PDBFixer: OK')"
python -c "import yaml; print('YAML: OK')"
python -c "import numpy; print('NumPy: OK')"
```

### Check Available Dependencies

Use PRISM's built-in dependency checker:

```python
import prism

# Check all dependencies
deps = prism.check_dependencies()
for name, available in deps.items():
    status = "✓" if available else "✗"
    print(f"{status} {name}")
```

Expected output:
```
✓ gromacs
✓ pdbfixer
✓ antechamber     # For GAFF/GAFF2
✗ openff          # Optional
✓ mdtraj          # Optional for analysis
✓ rdkit           # Optional for visualization
```

### List Available Force Fields

```bash
# Check installed protein force fields
prism --list-forcefields
```

### Test System Building

```bash
# Quick test with example files (if available)
prism test_protein.pdb test_ligand.mol2 -o test_output

# Check output structure
ls -lh test_output/
ls -lh test_output/GMX_PROLIG_MD/
```

## 🔧 Environment Variables

Set these environment variables for optimal performance:

```bash
# Add to ~/.bashrc or ~/.zshrc

# GROMACS path
export GMXRC=/opt/gromacs-2024.3/bin/GMXRC
source $GMXRC

# CUDA paths (adjust as needed)
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH

# PRISM path (optional, for direct usage)
export PRISM_HOME=/path/to/PRISM
export PATH=$PRISM_HOME/bin:$PATH
```

## 🐛 Troubleshooting

### Common Installation Issues

??? failure "GROMACS build fails with CUDA errors"

    **Solution**: Ensure CUDA paths are correct:
    ```bash
    # Find CUDA installation
    which nvcc
    ls /usr/local/cuda*
    
    # Update CMake command with correct paths
    cmake .. -DCUDA_TOOLKIT_ROOT_DIR=/correct/cuda/path ...
    ```

??? failure "ImportError: No module named 'prism'"

    **Solution**: Ensure you're in the correct conda environment:
    ```bash
    conda activate prism
    which python  # Should point to conda environment
    pip list | grep prism
    ```

??? failure "PDBFixer installation fails"

    **Solution**: Install from conda-forge channel:
    ```bash
    conda install -c conda-forge pdbfixer
    # If still fails, try:
    conda install -c conda-forge openmm pdbfixer
    ```

??? failure "ACPYPE not found"

    **Solution**: Install with pip in the correct environment:
    ```bash
    conda activate prism
    pip install acpype --upgrade
    ```

### Platform-Specific Issues

#### GPU Support on Linux

If GPU is not detected:

```bash
# Check NVIDIA driver
nvidia-smi

# Check CUDA installation
nvcc --version
ldconfig -p | grep cuda

# Reinstall CUDA toolkit if necessary
```

#### macOS Specific

For macOS users without CUDA:

```bash
# Build GROMACS without GPU support
cmake .. -DGMX_MPI=ON \
         -DGMX_BUILD_OWN_FFTW=ON \
         -DGMX_GPU=OFF \
         -DCMAKE_INSTALL_PREFIX=/opt/gromacs-2024.3
```

## 📚 Additional Resources

- [GROMACS Installation Guide](http://manual.gromacs.org/documentation/current/install-guide/index.html)
- [AmberTools Documentation](https://ambermd.org/AmberTools.php)
- [OpenFF Documentation](https://docs.openforcefield.org/)
- [PRISM GitHub Repository](https://github.com/your-username/PRISM)

## 🤝 Getting Help

If you encounter issues during installation:

1. Check the [Troubleshooting Guide](../user-guide/troubleshooting.md)
2. Search [GitHub Issues](https://github.com/AIB001/PRISM/issues)
3. Post on [GitHub Discussions](https://github.com/AIB001/PRISM/discussions)
4. Contact: [author email](mailto:zshi268@wisc.edu)

---

## ✨ Next Steps

After successful installation:

<div style="text-align: center; margin-top: 2rem;">
  <a href="quickstart/" class="md-button md-button--primary">Quick Start Tutorial</a>
  <a href="concepts/" class="md-button">Basic Concepts</a>
</div>