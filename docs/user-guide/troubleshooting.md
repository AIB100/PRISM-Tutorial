# Troubleshooting

Solutions for common errors when building systems, running simulations, or installing PRISM.

## Build Errors

### GROMACS Not Found

```
Error: GROMACS command not found
```

**Cause:** GROMACS is not installed, or its environment is not sourced in the current shell.

**Solution:**

```bash
# Source the GROMACS environment
source /path/to/gromacs/bin/GMXRC

# Verify it works
gmx --version

# Add to your shell profile for permanent access
echo 'source /path/to/gromacs/bin/GMXRC' >> ~/.bashrc
```

If GROMACS is installed but PRISM still cannot find it, specify the path explicitly:

```bash
prism protein.pdb ligand.mol2 -o output --gmx-command /opt/gromacs/bin/gmx
```

---

### Force Field Not Found

```
Error: Force field 'amber14sb' not found in GROMACS data directory
```

**Cause:** The protein force field you requested is not installed in your GROMACS data directory.

**Solution:**

```bash
# List force fields available in your GROMACS installation
prism --list-forcefields

# Install additional force fields (copies PRISM-bundled FFs to GROMACS)
prism --install-forcefields

# Or use a force field that GROMACS ships with
prism protein.pdb ligand.mol2 -o output --forcefield amber99sb
```

---

### Ligand Parameterization Failed

```
Error: Cannot assign atom types for ligand
```

**Cause:** The chosen force field generator could not assign parameters to your ligand. This usually happens with unusual functional groups or incorrect input geometry.

**Solution:**

1. **Check your ligand file** -- ensure it has correct bond orders, 3D coordinates, and explicit hydrogens.

2. **Try a different force field:**
    ```bash
    # OpenFF has broader chemical coverage
    prism protein.pdb ligand.sdf -o output --ligand-forcefield openff

    # Or try OPLS-AA via LigParGen
    prism protein.pdb ligand.mol2 -o output --ligand-forcefield opls
    ```

3. **Verify ligand format:** GAFF/GAFF2 require MOL2 with Tripos atom types. OpenFF works best with SDF.

4. **Check for missing hydrogens:**
    ```bash
    obabel ligand.mol2 -O ligand_h.mol2 -h  # Add hydrogens with OpenBabel
    ```

---

### Charge Assignment Failed

```
Error: antechamber failed to assign AM1-BCC charges
```

**Cause:** The semi-empirical charge calculation did not converge, often due to unusual chemistry, radicals, or very large ligands.

**Solution:**

```bash
# Specify the net charge explicitly
prism protein.pdb ligand.mol2 -o output --ligand-charge -1

# Or use Gaussian RESP charges for better accuracy
prism protein.pdb ligand.mol2 -o output --gaussian hf --nproc 8 --mem 4GB
```

If antechamber keeps failing, switch to OpenFF which uses a different charge method:

```bash
prism protein.pdb ligand.sdf -o output --ligand-forcefield openff
```

---

### Ion Addition Failed

```
Error: Cannot find solvent group for ion replacement
```

**Cause:** The system does not contain enough water molecules, or GROMACS cannot identify the solvent group automatically.

**Solution:**

- Increase the box size so more water is added:
    ```bash
    prism protein.pdb ligand.mol2 -o output --box-distance 2.0
    ```

- If the issue persists, check that solvation completed successfully by inspecting the `.gro` file for SOL residues.

---

### pdb2gmx Fails with Residue Errors

```
Fatal error: Residue 'XYZ' not found in residue topology database
```

**Cause:** The protein PDB contains non-standard residues (modified amino acids, ligands left in the file, or incorrect residue names).

**Solution:**

1. Remove non-protein atoms before running PRISM:
    ```bash
    grep "^ATOM" protein.pdb > protein_clean.pdb
    ```

2. For selenomethionine, rename MSE to MET:
    ```bash
    sed 's/MSE/MET/g; s/ SE / SD /' protein.pdb > protein_fixed.pdb
    ```

3. For other non-standard residues, use PDBFixer first:
    ```python
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    fixer = PDBFixer(filename='protein.pdb')
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    with open('protein_fixed.pdb', 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    ```

---

## Simulation Errors

### LINCS Warnings

```
LINCS WARNING: relative constraint deviation after LINCS
Step X  Bond between atoms X and X ...
```

**Cause:** A bond constraint was violated. Usually caused by bad initial contacts, too large a time step, or an unstable starting structure.

**Solution:**

1. **Run longer energy minimization:**
    ```bash
    prism protein.pdb ligand.mol2 -o output --em-steps 50000 --em-tolerance 100
    ```

2. **Use a smaller time step:**
    ```bash
    prism protein.pdb ligand.mol2 -o output --dt 0.001
    ```

3. **Increase box size** to reduce initial clashes:
    ```bash
    prism protein.pdb ligand.mol2 -o output --box-distance 2.0
    ```

4. **Check your ligand geometry** -- a ligand with distorted bond lengths will cause LINCS failures immediately. Re-optimize with RDKit or Avogadro before using.

---

### System Blew Up

```
Fatal error: Particle moved more than X nm
```

**Cause:** Atoms moved too far in a single step, causing the simulation to become numerically unstable.

**Solution:**

1. Re-run energy minimization with stricter tolerance
2. Extend equilibration time:
    ```bash
    prism protein.pdb ligand.mol2 -o output --nvt-ps 1000 --npt-ps 1000
    ```
3. Visualize the starting structure to check for overlapping atoms or ligand placed outside the protein
4. Reduce temperature coupling speed by editing `nvt.mdp`: set `tau_t = 1.0`

---

### GPU Errors

```
Cannot use GPU: no compatible devices detected
```

**Cause:** GROMACS was not compiled with GPU support, or CUDA drivers are missing.

**Solution:**

```bash
# Check GPU availability
nvidia-smi

# Check GROMACS GPU support
gmx mdrun -version | grep GPU

# If no GPU support, edit localrun.sh to remove GPU flags:
# Remove: -nb gpu -bonded gpu -pme gpu -gpu_id 0
```

For CPU-only runs, no changes are needed to PRISM -- just edit the generated `localrun.sh` script.

---

### Out of Memory

```
Cannot allocate memory / std::bad_alloc
```

**Cause:** The system is too large for available RAM, often during parameterization of very large ligands or very large solvated systems.

**Solution:**

- Reduce the box size: `--box-distance 1.0`
- Use dodecahedron box shape (fewer water molecules): `--box-shape dodecahedron`
- Reduce the number of OpenMP threads in `localrun.sh`: `-ntomp 4`
- On a cluster, request more memory: `#SBATCH --mem=64GB`

---

## Installation Errors

### ImportError: No module named 'prism'

**Cause:** PRISM is not installed in the active Python environment.

**Solution:**

```bash
# Ensure you are in the correct conda environment
conda activate prism

# Install PRISM in development mode
cd /path/to/PRISM
pip install -e .

# Verify
python -c "import prism; print('OK')"
```

---

### PDBFixer Installation Fails

**Cause:** PDBFixer depends on OpenMM, which requires conda-forge.

**Solution:**

```bash
conda install -c conda-forge pdbfixer
# If that fails:
conda install -c conda-forge openmm pdbfixer
```

---

### ACPYPE Not Found

**Cause:** ACPYPE is required for GAFF/GAFF2 but is not installed.

**Solution:**

```bash
pip install acpype
# Verify
acpype --help
```

---

### CUDA Version Mismatch

```
GROMACS was compiled with CUDA X but driver supports CUDA Y
```

**Cause:** The GROMACS binary was compiled against a different CUDA version than what your GPU driver provides.

**Solution:**

- Update your GPU driver to match the CUDA version GROMACS expects
- Or recompile GROMACS with the CUDA version your driver supports:
    ```bash
    cmake .. -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-12.4
    ```

---

## Advanced Errors

### PMF: Ligand Pulls Through Protein

**Cause:** The pulling direction passes through the protein body instead of pulling the ligand out of the binding pocket.

**Solution:**

Use the `--pull-vector` flag to specify two atom indices that define the correct unbinding direction (a protein atom near the pocket entrance and a ligand atom):

```bash
prism protein.pdb ligand.mol2 -o output --pmf --pull-vector 100 200
```

If you are unsure which atoms to use, PRISM's alignment module (used by default with `--pmf`) automatically selects an optimal pulling direction using Metropolis-Hastings sampling.

---

### PMF: WHAM Does Not Converge

**Cause:** Insufficient overlap between adjacent umbrella sampling windows.

**Solution:**

- Decrease the window spacing:
    ```bash
    prism ... --pmf --umbrella-spacing 0.08
    ```

- Increase the sampling time per window:
    ```bash
    prism ... --pmf --umbrella-time 20
    ```

- Check histogram overlap with `gmx wham -hist`

---

### MM/PBSA: gmx_MMPBSA Errors

**Cause:** Topology conversion or trajectory format issues.

**Solution:**

1. Ensure `gmx_MMPBSA` is installed:
    ```bash
    pip install gmx_MMPBSA
    ```

2. Use the `--gmx2amber` flag to handle topology conversion:
    ```bash
    prism protein.pdb ligand.mol2 -o output --mmpbsa --gmx2amber
    ```

3. If using a specific trajectory segment:
    ```bash
    prism protein.pdb ligand.mol2 -o output --mmpbsa --mmpbsa-traj 10
    ```

---

### REST2: Replica Exchange Rate Too Low

**Cause:** The temperature spacing between replicas is too large, causing poor exchange acceptance.

**Solution:**

- Increase the number of replicas:
    ```bash
    prism ... --rest2 --replica-number 24
    ```

- Narrow the temperature range:
    ```bash
    prism ... --rest2 --t-ref 310 --t-max 400
    ```

---

## General Tips

1. **Read the full error message.** PRISM prints detailed error messages with suggested fixes -- read them before searching online.

2. **Check file paths.** Most "file not found" errors are caused by running PRISM from the wrong directory.

3. **Try `--overwrite`.** If a previous run left partial output, use `-f` to force a clean rebuild:
    ```bash
    prism protein.pdb ligand.mol2 -o output -f
    ```

4. **Inspect intermediate files.** Look at `protein_clean.pdb` and the ligand `.itp` to catch issues early.

5. **Use `--export-config`** to see the exact parameters PRISM will use:
    ```bash
    prism --export-config debug_config.yaml
    ```

<div class="whats-next" markdown>

## What's Next

- [Return to the User Guide overview](index.md)
- [Review Output Files to understand what PRISM generates](output-files.md)
- [Open an issue on GitHub if your problem is not listed here](https://github.com/AIB001/PRISM/issues)

</div>
