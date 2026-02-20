# PMF API

Structure alignment and pulling direction optimization for PMF calculations.

!!! note "CLI vs Python API"
    The PMF workflow (steered MD, umbrella sampling, WHAM) is driven by the CLI with `--pmf`. The Python API documented here provides the alignment module used internally. For most users, the [PMF Calculations guide](../user-guide/pmf-calculations.md) covers everything needed.

## PMFAligner

Aligns a protein-ligand complex for PMF calculations by optimizing the pulling direction and rotating the system so the pull vector aligns with the Z-axis.

```python
from prism.pmf.alignment import PMFAligner

aligner = PMFAligner(pocket_cutoff=4.0, verbose=True)
```

### Parameters

- **pocket_cutoff** (float): Distance cutoff in Angstroms for defining pocket residues around the ligand. Default: `4.0`
- **verbose** (bool): Print detailed optimization output. Default: `True`

### Methods

#### align_for_pmf()

Align the protein-ligand complex for PMF pulling.

When `pullvec` is `None` (the default), the method uses Metropolis-Hastings optimization to find the pulling direction that maximizes clearance from binding pocket atoms. It then rotates the entire complex so this direction aligns with the +Z axis.

```python
results = aligner.align_for_pmf(
    protein_pdb="protein.pdb",
    ligand_file="ligand.mol2",
    output_dir="aligned_output",
    pullvec=None  # auto-detect optimal direction
)
```

**Parameters**:

- **protein_pdb** (str): Path to protein PDB file
- **ligand_file** (str): Path to ligand file (MOL2 or SDF)
- **output_dir** (str): Output directory for aligned structures
- **pullvec** (tuple of (int, int), optional): User-defined pull vector as `(protein_atom_index, ligand_atom_index)`. If `None`, uses MH-optimized direction.

**Returns**: dict with the following keys:

| Key | Type | Description |
| --- | --- | --- |
| `aligned_protein` | str | Path to aligned protein PDB |
| `aligned_ligand` | str | Path to aligned ligand file |
| `pull_vector` | ndarray | Normalized pull direction (unit vector) |
| `rotation_matrix` | ndarray | 3x3 rotation matrix applied |
| `rotation_center` | ndarray | Center of rotation (system geometric center) |
| `ligand_centroid` | ndarray | Ligand center of mass |
| `pocket_centroid` | ndarray | Pocket center of mass |
| `optimization` | dict | MH optimization details (auto mode only) |
| `pymol_script` | str | Path to PyMOL visualization script (auto mode only) |

---

## How It Works

The alignment process:

1. **Parse structures** - Read protein PDB and ligand MOL2/SDF
2. **Identify pocket** - Find protein heavy atoms within `pocket_cutoff` of any ligand atom (vectorized distance calculation)
3. **Optimize direction** - Run Metropolis-Hastings simulated annealing to find the pulling direction with maximum clearance from pocket atoms
4. **Rotate complex** - Apply rotation matrix to align the optimized pull vector with +Z
5. **Shift coordinates** - Ensure all coordinates are positive (required by GROMACS)
6. **Write output** - Save aligned protein and ligand, plus a PyMOL visualization script

---

## Examples

### Auto-Optimized Direction

```python
from prism.pmf.alignment import PMFAligner

aligner = PMFAligner(pocket_cutoff=4.0, verbose=True)

results = aligner.align_for_pmf(
    "protein.pdb",
    "ligand.mol2",
    "pmf_aligned"
)

print(f"Aligned protein: {results['aligned_protein']}")
print(f"Aligned ligand: {results['aligned_ligand']}")
print(f"Pull vector: {results['pull_vector']}")

# Check optimization quality
opt = results['optimization']
print(f"Pocket residues: {opt['pocket_residues']}")
print(f"Energy improvement: {opt['improvement_pct']:.1f}%")
```

### User-Defined Pull Vector

```python
aligner = PMFAligner(verbose=True)

results = aligner.align_for_pmf(
    "protein.pdb",
    "ligand.mol2",
    "pmf_aligned",
    pullvec=(100, 200)  # protein atom 100, ligand atom 200
)
```

### Integration with Builder (Internal)

The `PRISMBuilder` uses `PMFAligner` internally when the `--pmf` flag is set:

```python
# This is what happens inside builder.py when --pmf is used:
from prism.pmf.alignment import PMFAligner

aligner = PMFAligner(pocket_cutoff=4.0, verbose=True)
alignment_results = aligner.align_for_pmf(
    protein_pdb, ligand_file, output_dir, pullvec
)

# Builder then uses aligned structures for system building
aligned_protein = alignment_results['aligned_protein']
aligned_ligand = alignment_results['aligned_ligand']
```

---

## CLI Usage

For most workflows, use the CLI flags instead of the Python API:

```bash
# Auto-optimized pulling direction (recommended)
prism protein.pdb ligand.mol2 -o output --pmf

# User-defined pulling direction
prism protein.pdb ligand.mol2 -o output --pmf --pull-vector 100 200
```

See the [PMF Calculations guide](../user-guide/pmf-calculations.md) for the complete CLI workflow including umbrella sampling and WHAM analysis.

---

## Related

- [User Guide: PMF Calculations](../user-guide/pmf-calculations.md)
- [Tutorial: PMF Calculations](../tutorials/pmf-tutorial.md)
