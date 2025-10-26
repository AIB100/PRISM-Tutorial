# Batch Processing Tutorial

## Overview

Process multiple protein-ligand complexes efficiently using PRISM's Python API and automation features. This tutorial demonstrates high-throughput virtual screening workflows.

**Time:** 60-90 minutes
**Level:** Intermediate

## Prerequisites

- Python programming basics
- Completed [Basic Tutorial](basic-tutorial.md)
- Multiple ligand structures (MOL2 or SDF format)
- Understanding of parallel computing concepts

## Use Cases

1. **Virtual Screening**: Build and simulate 100+ protein-ligand complexes
2. **Force Field Comparison**: Test same system with multiple force fields
3. **Mutation Studies**: Compare wild-type vs mutant proteins
4. **Trajectory Processing**: Batch analyze multiple MD trajectories

---

## Workflow 1: Virtual Screening Pipeline

### Step 1: Organize Input Files

```bash
mkdir batch_screening
cd batch_screening

# Create directory structure
mkdir -p ligands results logs
```

Expected directory structure:
```
batch_screening/
├── protein.pdb              # Target protein
├── ligands/                 # Ligand library
│   ├── compound_001.mol2
│   ├── compound_002.mol2
│   ├── compound_003.mol2
│   └── ...
├── results/                 # Output will go here
└── logs/                    # Build logs
```

### Step 2: Batch System Building

Create `batch_build.py`:

```python
#!/usr/bin/env python3
import prism
from pathlib import Path
import multiprocessing as mp
from datetime import datetime
import json

# Configuration
PROTEIN = "protein.pdb"
LIGAND_DIR = Path("ligands")
OUTPUT_DIR = Path("results")
LOG_DIR = Path("logs")
N_PROCESSES = 4  # Adjust based on your CPU cores

# Create directories
OUTPUT_DIR.mkdir(exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

def build_single_system(ligand_path):
    """Build a single protein-ligand system"""
    ligand_name = ligand_path.stem
    output_path = OUTPUT_DIR / ligand_name
    log_file = LOG_DIR / f"{ligand_name}.log"

    start_time = datetime.now()

    try:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Building {ligand_name}...")

        # Build system
        system = prism.system(
            PROTEIN,
            str(ligand_path),
            output_dir=str(output_path),
            ligand_forcefield="gaff",
            forcefield="amber99sb-ildn",
            water_model="tip3p",
            production_ns=100  # Shorter for screening
        )
        system.build()

        # Record success
        elapsed = (datetime.now() - start_time).total_seconds()
        result = {
            "ligand": ligand_name,
            "status": "SUCCESS",
            "time_seconds": elapsed,
            "output_dir": str(output_path)
        }

        # Save log
        with open(log_file, 'w') as f:
            f.write(f"Build successful in {elapsed:.1f}s\n")
            f.write(f"Output: {output_path}\n")

        return result

    except Exception as e:
        elapsed = (datetime.now() - start_time).total_seconds()
        result = {
            "ligand": ligand_name,
            "status": "FAILED",
            "time_seconds": elapsed,
            "error": str(e)
        }

        # Save error log
        with open(log_file, 'w') as f:
            f.write(f"Build failed after {elapsed:.1f}s\n")
            f.write(f"Error: {e}\n")

        return result

def main():
    # Find all ligand files
    ligand_files = sorted(LIGAND_DIR.glob("*.mol2")) + sorted(LIGAND_DIR.glob("*.sdf"))

    print(f"\n{'='*60}")
    print(f"PRISM Batch Processing")
    print(f"{'='*60}")
    print(f"Protein: {PROTEIN}")
    print(f"Ligands found: {len(ligand_files)}")
    print(f"Parallel processes: {N_PROCESSES}")
    print(f"{'='*60}\n")

    if len(ligand_files) == 0:
        print("ERROR: No ligand files found!")
        return

    # Process in parallel
    start_time = datetime.now()

    with mp.Pool(processes=N_PROCESSES) as pool:
        results = pool.map(build_single_system, ligand_files)

    total_time = (datetime.now() - start_time).total_seconds()

    # Summary statistics
    successes = [r for r in results if r['status'] == 'SUCCESS']
    failures = [r for r in results if r['status'] == 'FAILED']

    print(f"\n{'='*60}")
    print(f"BATCH BUILD SUMMARY")
    print(f"{'='*60}")
    print(f"Total ligands: {len(results)}")
    print(f"Successful: {len(successes)} ({100*len(successes)/len(results):.1f}%)")
    print(f"Failed: {len(failures)} ({100*len(failures)/len(results):.1f}%)")
    print(f"Total time: {total_time/60:.1f} minutes")
    print(f"Average time per system: {total_time/len(results):.1f} seconds")
    print(f"{'='*60}\n")

    # Save detailed results
    with open('build_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Print failures if any
    if failures:
        print("Failed ligands:")
        for f in failures:
            print(f"  - {f['ligand']}: {f['error']}")

    print(f"\nDetailed results saved to: build_results.json")

if __name__ == "__main__":
    main()
```

Run the batch building:

```bash
python batch_build.py
```

### Step 3: Monitor Progress

While building is running, monitor in another terminal:

```bash
# Watch progress
watch -n 5 'ls results/ | wc -l'

# Check latest logs
tail -f logs/*.log

# Check for errors
grep -i error logs/*.log
```

### Step 4: Batch Simulation Submission

#### Option A: SLURM Cluster

Create `submit_all_slurm.sh`:

```bash
#!/bin/bash
# Submit all systems to SLURM cluster

for system_dir in results/*/GMX_PROLIG_MD; do
    ligand_name=$(basename $(dirname $system_dir))

    # Create individual SLURM script
    cat > ${system_dir}/submit.slurm << EOF
#!/bin/bash
#SBATCH --job-name=${ligand_name}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --partition=gpu
#SBATCH --output=${ligand_name}_%j.out
#SBATCH --error=${ligand_name}_%j.err

# Load modules
module load gromacs/2024.3

# Run simulation
cd ${system_dir}
bash localrun.sh

echo "Simulation completed for ${ligand_name}"
EOF

    # Submit job
    cd $system_dir
    sbatch submit.slurm
    cd - > /dev/null

    echo "Submitted: $ligand_name"
done

echo "All jobs submitted!"
```

Run:
```bash
chmod +x submit_all_slurm.sh
./submit_all_slurm.sh
```

#### Option B: Local Workstation (Sequential)

```bash
#!/bin/bash
# run_all_local.sh

for system_dir in results/*/GMX_PROLIG_MD; do
    ligand_name=$(basename $(dirname $system_dir))
    echo "Running simulation for: $ligand_name"

    cd $system_dir
    bash localrun.sh > ${ligand_name}_md.log 2>&1
    cd - > /dev/null

    echo "Completed: $ligand_name"
done
```

### Step 5: Batch Trajectory Analysis

After simulations complete, analyze all trajectories:

```python
#!/usr/bin/env python3
# batch_analyze.py

import prism
from pathlib import Path
import pandas as pd
import multiprocessing as mp

def analyze_trajectory(traj_dir):
    """Analyze a single trajectory"""
    ligand_name = traj_dir.parent.parent.name

    try:
        # Paths
        topology = str(traj_dir / "md.tpr")
        trajectory = str(traj_dir / "md.xtc")
        output_dir = str(traj_dir / "analysis")

        print(f"Analyzing {ligand_name}...")

        # Run analysis
        analyzer = prism.analyze_trajectory(
            topology=topology,
            trajectory=trajectory,
            ligand_resname="LIG",
            output_dir=output_dir
        )

        # Extract key metrics
        results = {
            "ligand": ligand_name,
            "rmsd_avg": analyzer.results.get('rmsd_avg', None),
            "rmsd_std": analyzer.results.get('rmsd_std', None),
            "contacts_avg": analyzer.results.get('contacts_avg', None),
            "hbonds_avg": analyzer.results.get('hbonds_avg', None),
            "status": "SUCCESS"
        }

        return results

    except Exception as e:
        return {
            "ligand": ligand_name,
            "status": "FAILED",
            "error": str(e)
        }

def main():
    # Find all trajectory directories
    traj_dirs = list(Path("results").glob("*/GMX_PROLIG_MD/prod"))
    traj_dirs = [d for d in traj_dirs if (d / "md.xtc").exists()]

    print(f"Found {len(traj_dirs)} trajectories to analyze")

    # Analyze in parallel
    with mp.Pool(processes=4) as pool:
        results = pool.map(analyze_trajectory, traj_dirs)

    # Create summary DataFrame
    df = pd.DataFrame(results)
    df.to_csv("analysis_summary.csv", index=False)

    print("\nAnalysis Summary:")
    print(df)
    print(f"\nResults saved to: analysis_summary.csv")

if __name__ == "__main__":
    main()
```

Run analysis:
```bash
python batch_analyze.py
```

---

## Workflow 2: Force Field Comparison

Compare same system with multiple force fields:

```python
#!/usr/bin/env python3
# compare_forcefields.py

import prism

# Protein and ligand
protein = "protein.pdb"
ligand = "ligand.mol2"

# Force fields to compare
force_fields = [
    {"name": "GAFF", "ligand_ff": "gaff"},
    {"name": "GAFF2", "ligand_ff": "gaff2"},
    {"name": "OpenFF", "ligand_ff": "openff"},
    {"name": "OPLS-AA", "ligand_ff": "opls"}
]

# Build systems
for ff in force_fields:
    print(f"\nBuilding with {ff['name']}...")

    output_dir = f"comparison_{ff['name'].lower()}"

    try:
        system = prism.system(
            protein, ligand,
            output_dir=output_dir,
            ligand_forcefield=ff['ligand_ff'],
            forcefield="amber99sb-ildn"
        )
        system.build()
        print(f"  SUCCESS: {output_dir}")

    except Exception as e:
        print(f"  FAILED: {e}")

print("\nAll force field systems built!")
print("Next: Run MD simulations and compare results")
```

---

## Workflow 3: Batch Trajectory Processing

Process multiple DCD trajectories:

```python
#!/usr/bin/env python3
# batch_process_trajectories.py

import prism
from pathlib import Path

# Configuration
INPUT_DIR = Path("raw_trajectories")
OUTPUT_DIR = Path("processed_trajectories")
TOPOLOGY = "system.tpr"

# Find all trajectories
trajectories = list(INPUT_DIR.glob("*.dcd"))
print(f"Found {len(trajectories)} DCD files")

# Process in batch
processed = prism.batch_process_trajectories(
    input_trajectories=[str(t) for t in trajectories],
    output_dir=str(OUTPUT_DIR),
    topology_file=TOPOLOGY,
    prefix="processed_",
    center_selection="Protein",
    output_selection="System",
    pbc_method="mol"
)

print(f"\nProcessed {len(processed)} trajectories:")
for p in processed:
    print(f"  - {Path(p).name}")
```

---

## Advanced: Custom Analysis Pipeline

```python
#!/usr/bin/env python3
# custom_pipeline.py

import prism
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class ScreeningAnalyzer:
    """Custom analyzer for virtual screening results"""

    def __init__(self, results_dir):
        self.results_dir = Path(results_dir)
        self.data = []

    def analyze_all(self):
        """Analyze all completed simulations"""
        systems = list(self.results_dir.glob("*/GMX_PROLIG_MD"))

        for system_dir in systems:
            ligand_name = system_dir.parent.name
            self.analyze_system(ligand_name, system_dir)

    def analyze_system(self, name, system_dir):
        """Analyze a single system"""
        traj_file = system_dir / "prod" / "md.xtc"
        topo_file = system_dir / "prod" / "md.tpr"

        if not traj_file.exists():
            return

        try:
            # Load and analyze
            analyzer = prism.TrajAnalysis(
                str(topo_file),
                str(traj_file),
                ligand_resname="LIG"
            )

            # Extract metrics
            metrics = {
                "ligand": name,
                "rmsd_protein": analyzer.calc_rmsd(selection="protein"),
                "rmsd_ligand": analyzer.calc_rmsd(selection="resname LIG"),
                "contacts": analyzer.calc_contacts(),
                "hbonds": analyzer.calc_hbonds(),
                "sasa_ligand": analyzer.calc_sasa(selection="resname LIG")
            }

            self.data.append(metrics)

        except Exception as e:
            print(f"Failed to analyze {name}: {e}")

    def rank_ligands(self):
        """Rank ligands by stability"""
        df = pd.DataFrame(self.data)

        # Calculate stability score (lower RMSD + more contacts = better)
        df['stability_score'] = (
            -df['rmsd_ligand'].mean(axis=1) +
            df['contacts'].mean(axis=1) * 0.1 +
            df['hbonds'].mean(axis=1) * 0.5
        )

        # Rank
        df_ranked = df.sort_values('stability_score', ascending=False)
        return df_ranked

    def plot_summary(self, output="screening_summary.png"):
        """Create summary plot"""
        df_ranked = self.rank_ligands()

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # RMSD distribution
        axes[0, 0].boxplot([df_ranked['rmsd_ligand']])
        axes[0, 0].set_title("Ligand RMSD Distribution")
        axes[0, 0].set_ylabel("RMSD (nm)")

        # Contacts
        axes[0, 1].scatter(df_ranked['contacts'].mean(axis=1),
                           df_ranked['stability_score'])
        axes[0, 1].set_xlabel("Average Contacts")
        axes[0, 1].set_ylabel("Stability Score")
        axes[0, 1].set_title("Contacts vs Stability")

        # H-bonds
        axes[1, 0].hist(df_ranked['hbonds'].mean(axis=1), bins=20)
        axes[1, 0].set_xlabel("Average H-bonds")
        axes[1, 0].set_ylabel("Frequency")
        axes[1, 0].set_title("H-bond Distribution")

        # Top candidates
        top_10 = df_ranked.head(10)
        axes[1, 1].barh(range(len(top_10)), top_10['stability_score'])
        axes[1, 1].set_yticks(range(len(top_10)))
        axes[1, 1].set_yticklabels(top_10['ligand'])
        axes[1, 1].set_xlabel("Stability Score")
        axes[1, 1].set_title("Top 10 Candidates")

        plt.tight_layout()
        plt.savefig(output, dpi=300)
        print(f"Summary plot saved to: {output}")

# Usage
analyzer = ScreeningAnalyzer("results")
analyzer.analyze_all()
df_ranked = analyzer.rank_ligands()
df_ranked.to_csv("ranked_ligands.csv")
analyzer.plot_summary()

print("\nTop 10 ligands:")
print(df_ranked.head(10)[['ligand', 'stability_score']])
```

---

## Best Practices

### 1. Resource Management

```python
# Don't run too many parallel builds
# Rule of thumb: N_PROCESSES = CPU_cores / 2
N_PROCESSES = max(1, mp.cpu_count() // 2)
```

### 2. Error Handling

Always use try-except blocks and log errors:

```python
def safe_build(ligand):
    try:
        return build_system(ligand)
    except Exception as e:
        log_error(ligand, e)
        return None
```

### 3. Checkpointing

Save progress periodically:

```python
# Save results after each ligand
with open(f'checkpoint_{ligand_name}.json', 'w') as f:
    json.dump(result, f)
```

### 4. Disk Space Monitoring

```bash
# Check disk usage
du -sh results/
df -h .
```

---

## Troubleshooting

### Memory Issues

Reduce parallel processes:
```python
N_PROCESSES = 2  # Instead of 4
```

### Failed Builds

Resume from failures:
```python
# Load previous results
with open('build_results.json') as f:
    previous = json.load(f)

# Find failures
failed = [r['ligand'] for r in previous if r['status'] == 'FAILED']

# Retry failures
for ligand_name in failed:
    retry_build(ligand_name)
```

### Slow I/O

Use faster storage (SSD) for output:
```python
OUTPUT_DIR = Path("/scratch/fast_storage/results")
```

---

## Next Steps

- Integrate with docking workflows
- Add free energy calculations (PMF)
- Implement machine learning for hit prediction
- Scale to HPC clusters for thousands of ligands

---

## Related Resources

- [Basic Tutorial](basic-tutorial.md)
- [PMF Tutorial](pmf-tutorial.md)
- [API Documentation](../api/index.md)
