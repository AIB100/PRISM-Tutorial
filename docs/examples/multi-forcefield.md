# Multiple Force Fields Comparison Example

Compare protein-ligand binding across different force fields to assess force field dependence and validate results.

## Overview

This example demonstrates:

- Building the same system with multiple force fields
- Running parallel simulations
- Comparing structural and energetic properties
- Statistical analysis of force field effects

**Complexity**: Intermediate  
**Time**: 2-3 hours + simulation time  
**Prerequisites**: Completed [Simple Example](simple.md)

---

## Supported Force Field Combinations

### Ligand Force Fields

| Force Field | Description | Best For | Dependencies |
|------------|-------------|----------|--------------|
| **GAFF** | General AMBER Force Field | General purpose | AmberTools, ACPYPE |
| **GAFF2** | Improved GAFF | Better for drug-like | AmberTools, ACPYPE |
| **OpenFF** | Open Force Field | Modern, data-driven | OpenFF Toolkit |
| **OPLS-AA** | OPLS All-Atom | Proteins & small molecules | LigParGen (web) |
| **CGenFF** | CHARMM General FF | CHARMM compatibility | CGenFF program |
| **MMFF** | Merck Molecular FF | Quick screening | SwissParam (web) |

### Protein Force Fields

- amber99sb (default)
- amber99sb-ildn (improved side chains)
- amber14sb (recommended)
- charmm36 (for membranes)
- oplsaa (consistent with OPLS ligands)

---

## Complete Example Script

```python
#!/usr/bin/env python3
"""
Multi-Force Field Comparison Example
Compare T4 Lysozyme-Benzene with different force fields
"""

import prism
import os
import multiprocessing as mp
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

class ForceFieldComparator:
    """Compare multiple force field combinations"""
    
    def __init__(self, protein, ligand, base_dir="ff_comparison"):
        self.protein = protein
        self.ligand = ligand
        self.base_dir = Path(base_dir)
        self.base_dir.mkdir(exist_ok=True)
        self.results = []
        
    def build_all_systems(self, force_fields):
        """Build systems with different force fields"""
        print("="*70)
        print("Building Systems with Multiple Force Fields")
        print("="*70)
        
        for i, ff_config in enumerate(force_fields, 1):
            print(f"\n[{i}/{len(force_fields)}] {ff_config['name']}")
            print("-" * 70)
            
            output_dir = self.base_dir / ff_config['dir_name']
            
            try:
                system = prism.system(
                    self.protein,
                    self.ligand,
                    output_dir=str(output_dir),
                    ligand_forcefield=ff_config.get('ligand_ff', 'gaff'),
                    forcefield=ff_config.get('protein_ff', 'amber99sb'),
                    water_model=ff_config.get('water', 'tip3p'),
                    production_ns=100  # Shorter for comparison
                )
                
                system.build()
                
                result = {
                    'name': ff_config['name'],
                    'dir': str(output_dir),
                    'ligand_ff': ff_config.get('ligand_ff', 'gaff'),
                    'protein_ff': ff_config.get('protein_ff', 'amber99sb'),
                    'status': 'SUCCESS',
                    'output_dir': output_dir
                }
                
                print(f"✓ SUCCESS: {output_dir}")
                
            except Exception as e:
                result = {
                    'name': ff_config['name'],
                    'status': 'FAILED',
                    'error': str(e)
                }
                print(f"✗ FAILED: {e}")
            
            self.results.append(result)
        
        # Summary
        successes = [r for r in self.results if r['status'] == 'SUCCESS']
        print(f"\n{'='*70}")
        print(f"Build Summary: {len(successes)}/{len(self.results)} succeeded")
        print(f"{'='*70}\n")
        
        return self.results
    
    def analyze_all_trajectories(self):
        """Analyze all completed simulations"""
        print("\nAnalyzing Trajectories...")
        
        analysis_results = []
        
        for result in self.results:
            if result['status'] != 'SUCCESS':
                continue
                
            traj_dir = Path(result['output_dir']) / "GMX_PROLIG_MD" / "prod"
            traj_file = traj_dir / "md.xtc"
            topo_file = traj_dir / "md.tpr"
            
            if not traj_file.exists():
                print(f"  Skipping {result['name']}: trajectory not found")
                continue
            
            print(f"  Analyzing {result['name']}...")
            
            try:
                # Basic RMSD analysis with MDTraj
                import mdtraj as md
                
                traj = md.load(str(traj_file), top=str(topo_file))
                
                # Calculate RMSDs
                protein_atoms = traj.topology.select("protein and name CA")
                ligand_atoms = traj.topology.select("resname LIG")
                
                rmsd_protein = md.rmsd(traj, traj, 0, atom_indices=protein_atoms) * 10  # nm to Å
                rmsd_ligand = md.rmsd(traj, traj, 0, atom_indices=ligand_atoms) * 10
                
                analysis_results.append({
                    'name': result['name'],
                    'ligand_ff': result['ligand_ff'],
                    'protein_ff': result['protein_ff'],
                    'rmsd_protein_mean': rmsd_protein.mean(),
                    'rmsd_protein_std': rmsd_protein.std(),
                    'rmsd_ligand_mean': rmsd_ligand.mean(),
                    'rmsd_ligand_std': rmsd_ligand.std(),
                    'rmsd_protein': rmsd_protein,
                    'rmsd_ligand': rmsd_ligand
                })
                
            except Exception as e:
                print(f"  Failed to analyze {result['name']}: {e}")
        
        return analysis_results
    
    def plot_comparison(self, analysis_results, output="ff_comparison.png"):
        """Create comparison plots"""
        if not analysis_results:
            print("No analysis results to plot")
            return
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        names = [r['name'] for r in analysis_results]
        
        # Plot 1: Protein RMSD mean
        protein_rmsd_means = [r['rmsd_protein_mean'] for r in analysis_results]
        protein_rmsd_stds = [r['rmsd_protein_std'] for r in analysis_results]
        
        axes[0, 0].bar(range(len(names)), protein_rmsd_means, yerr=protein_rmsd_stds)
        axes[0, 0].set_xticks(range(len(names)))
        axes[0, 0].set_xticklabels(names, rotation=45, ha='right')
        axes[0, 0].set_ylabel('RMSD (Å)')
        axes[0, 0].set_title('Protein Backbone RMSD')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Ligand RMSD mean
        ligand_rmsd_means = [r['rmsd_ligand_mean'] for r in analysis_results]
        ligand_rmsd_stds = [r['rmsd_ligand_std'] for r in analysis_results]
        
        axes[0, 1].bar(range(len(names)), ligand_rmsd_means, yerr=ligand_rmsd_stds)
        axes[0, 1].set_xticks(range(len(names)))
        axes[0, 1].set_xticklabels(names, rotation=45, ha='right')
        axes[0, 1].set_ylabel('RMSD (Å)')
        axes[0, 1].set_title('Ligand RMSD')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Protein RMSD time series
        for r in analysis_results:
            time = np.arange(len(r['rmsd_protein'])) * 0.5  # Assuming 500 ps output
            axes[1, 0].plot(time, r['rmsd_protein'], label=r['name'], alpha=0.7)
        
        axes[1, 0].set_xlabel('Time (ns)')
        axes[1, 0].set_ylabel('Protein RMSD (Å)')
        axes[1, 0].set_title('Protein RMSD Over Time')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Ligand RMSD time series
        for r in analysis_results:
            time = np.arange(len(r['rmsd_ligand'])) * 0.5
            axes[1, 1].plot(time, r['rmsd_ligand'], label=r['name'], alpha=0.7)
        
        axes[1, 1].set_xlabel('Time (ns)')
        axes[1, 1].set_ylabel('Ligand RMSD (Å)')
        axes[1, 1].set_title('Ligand RMSD Over Time')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"\nComparison plot saved: {output}")
    
    def generate_report(self, analysis_results):
        """Generate comparison report"""
        df = pd.DataFrame([{
            'Force Field': r['name'],
            'Ligand FF': r['ligand_ff'],
            'Protein FF': r['protein_ff'],
            'Protein RMSD (Å)': f"{r['rmsd_protein_mean']:.2f} ± {r['rmsd_protein_std']:.2f}",
            'Ligand RMSD (Å)': f"{r['rmsd_ligand_mean']:.2f} ± {r['rmsd_ligand_std']:.2f}"
        } for r in analysis_results])
        
        report_file = self.base_dir / "comparison_report.csv"
        df.to_csv(report_file, index=False)
        
        print("\n" + "="*70)
        print("Force Field Comparison Report")
        print("="*70)
        print(df.to_string(index=False))
        print(f"\nReport saved: {report_file}")

def main():
    # Input files
    protein = "protein.pdb"
    ligand = "ligand.mol2"
    
    # Force field configurations to test
    force_fields = [
        {
            'name': 'GAFF + AMBER99SB',
            'dir_name': 'gaff_amber99sb',
            'ligand_ff': 'gaff',
            'protein_ff': 'amber99sb',
            'water': 'tip3p'
        },
        {
            'name': 'GAFF2 + AMBER99SB',
            'dir_name': 'gaff2_amber99sb',
            'ligand_ff': 'gaff2',
            'protein_ff': 'amber99sb',
            'water': 'tip3p'
        },
        {
            'name': 'OpenFF + AMBER14SB',
            'dir_name': 'openff_amber14sb',
            'ligand_ff': 'openff',
            'protein_ff': 'amber14sb',
            'water': 'tip3p'
        },
        {
            'name': 'GAFF + AMBER14SB',
            'dir_name': 'gaff_amber14sb',
            'ligand_ff': 'gaff',
            'protein_ff': 'amber14sb',
            'water': 'tip3p'
        }
    ]
    
    # Create comparator
    comparator = ForceFieldComparator(protein, ligand)
    
    # Build all systems
    build_results = comparator.build_all_systems(force_fields)
    
    print("\n" + "="*70)
    print("Next Steps:")
    print("="*70)
    print("Run simulations for each force field:")
    for result in build_results:
        if result['status'] == 'SUCCESS':
            sim_dir = Path(result['output_dir']) / "GMX_PROLIG_MD"
            print(f"  cd {sim_dir} && bash localrun.sh &")
    
    print("\nAfter all simulations complete, run analysis:")
    print("  python run_analysis.py")

if __name__ == "__main__":
    main()
```

Create separate analysis script `run_analysis.py`:

```python
#!/usr/bin/env python3
"""Analyze force field comparison results"""

from pathlib import Path
import sys

# Load the comparator with existing results
sys.path.insert(0, '.')
from run_multi_forcefield import ForceFieldComparator

def main():
    comparator = ForceFieldComparator("protein.pdb", "ligand.mol2")
    
    # Load build results
    comparator.results = [
        {'name': 'GAFF + AMBER99SB', 'status': 'SUCCESS', 
         'output_dir': Path('ff_comparison/gaff_amber99sb'),
         'ligand_ff': 'gaff', 'protein_ff': 'amber99sb'},
        {'name': 'GAFF2 + AMBER99SB', 'status': 'SUCCESS',
         'output_dir': Path('ff_comparison/gaff2_amber99sb'),
         'ligand_ff': 'gaff2', 'protein_ff': 'amber99sb'},
        {'name': 'OpenFF + AMBER14SB', 'status': 'SUCCESS',
         'output_dir': Path('ff_comparison/openff_amber14sb'),
         'ligand_ff': 'openff', 'protein_ff': 'amber14sb'},
        {'name': 'GAFF + AMBER14SB', 'status': 'SUCCESS',
         'output_dir': Path('ff_comparison/gaff_amber14sb'),
         'ligand_ff': 'gaff', 'protein_ff': 'amber14sb'}
    ]
    
    # Analyze all trajectories
    analysis_results = comparator.analyze_all_trajectories()
    
    # Plot comparison
    comparator.plot_comparison(analysis_results)
    
    # Generate report
    comparator.generate_report(analysis_results)

if __name__ == "__main__":
    main()
```

---

## Running the Example

### Step 1: Build All Systems

```bash
python run_multi_forcefield.py
```

This will create:
```
ff_comparison/
├── gaff_amber99sb/
├── gaff2_amber99sb/
├── openff_amber14sb/
└── gaff_amber14sb/
```

### Step 2: Run Simulations

Option A - Sequential:
```bash
cd ff_comparison/gaff_amber99sb/GMX_PROLIG_MD && bash localrun.sh
cd ../../../ff_comparison/gaff2_amber99sb/GMX_PROLIG_MD && bash localrun.sh
# etc...
```

Option B - Parallel (if you have multiple GPUs):
```bash
cd ff_comparison/gaff_amber99sb/GMX_PROLIG_MD && bash localrun.sh &
cd ../../gaff2_amber99sb/GMX_PROLIG_MD && bash localrun.sh &
cd ../../openff_amber14sb/GMX_PROLIG_MD && bash localrun.sh &
cd ../../gaff_amber14sb/GMX_PROLIG_MD && bash localrun.sh &
wait
```

### Step 3: Analyze and Compare

After all simulations complete:

```bash
python run_analysis.py
```

---

## Expected Results

### Typical RMSD Ranges

For T4 Lysozyme-Benzene system:

| Force Field | Protein RMSD | Ligand RMSD |
|------------|--------------|-------------|
| GAFF + AMBER99SB | 1.2 ± 0.3 Å | 1.5 ± 0.8 Å |
| GAFF2 + AMBER99SB | 1.1 ± 0.3 Å | 1.3 ± 0.7 Å |
| OpenFF + AMBER14SB | 1.0 ± 0.2 Å | 1.4 ± 0.6 Å |
| GAFF + AMBER14SB | 1.1 ± 0.2 Å | 1.6 ± 0.9 Å |

### Interpretation

- **Protein RMSD < 2 Å**: Good stability, force field performs well
- **Ligand RMSD < 2 Å**: Ligand stays bound, good sampling
- **Small differences between FFs**: Results are robust
- **Large differences (>1 Å)**: Force field dependent, need validation

---

## Advanced Analysis

### Statistical Significance

```python
from scipy import stats

# Compare GAFF vs OpenFF
gaff_rmsd = analysis_results[0]['rmsd_ligand']
openff_rmsd = analysis_results[2]['rmsd_ligand']

# T-test
t_stat, p_value = stats.ttest_ind(gaff_rmsd, openff_rmsd)
print(f"P-value: {p_value:.4f}")

if p_value < 0.05:
    print("Significant difference between force fields")
else:
    print("No significant difference")
```

---

## Best Practices

1. **Use Same Settings**: Keep all parameters except force field constant
2. **Multiple Replicates**: Run 3+ simulations per force field
3. **Sufficient Sampling**: At least 100 ns per simulation
4. **Statistical Testing**: Use t-tests or ANOVA for comparison
5. **Validate**: Compare with experimental data when available

---

## Next Steps

- [Custom Configuration Example](custom-config.md)
- [PMF Calculations](../tutorials/pmf-tutorial.md)
- [Batch Processing](../tutorials/batch-tutorial.md)
