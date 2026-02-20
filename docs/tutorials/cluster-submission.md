# Cluster Submission Templates

Ready-to-use Slurm job scripts for running PRISM workflows on HPC clusters. All scripts include GPU support, checkpoint resume, and diagnostic logging.

!!! tip "Copy or Download"
    Click the copy button (:material-content-copy:) at the top-right corner of any code block to copy the script to your clipboard.

---

## Generic Template

A minimal template that works with any PRISM workflow. Replace the placeholder values marked with `<...>` before submitting.

```bash title="submit.sh"
#!/bin/bash
#SBATCH -J PRISM
#SBATCH -p <your_partition>           # e.g., gpu, gpu_a100, normal
#SBATCH --time=168:00:00              # 7 days walltime
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -o slurm-%j.out
#SBATCH -e slurm-%j.err

######################################################
# ENVIRONMENT SETUP
######################################################

# Option 1: Module system (most clusters)
module purge
module load gromacs/2024.3          # adjust version to your cluster

# Option 2: Manual environment (if no module system)
# export LD_LIBRARY_PATH=/path/to/libs:$LD_LIBRARY_PATH
# source /path/to/GMXRC

######################################################
# DIAGNOSTICS
######################################################

echo "=========================================="
echo "Start time:             $(date)"
echo "Job ID:                 $SLURM_JOB_ID"
echo "Node:                   $(hostname)"
echo "Working directory:      $(pwd)"
echo "CUDA_VISIBLE_DEVICES:   $CUDA_VISIBLE_DEVICES"
echo "GROMACS version:        $(gmx --version 2>&1 | head -1)"
echo "=========================================="

######################################################
# SIMULATION
######################################################

cd $SLURM_SUBMIT_DIR/GMX_PROLIG_MD   # adjust path to your output directory

bash localrun.sh

echo "End time: $(date)"
```

Submit with:

```bash
sbatch submit.sh
```

---

## Standard MD

For a standard protein-ligand MD simulation built with PRISM:

```bash title="submit_md.sh"
#!/bin/bash
#SBATCH -J PRISM_MD
#SBATCH -p <your_partition>
#SBATCH --time=72:00:00               # 3 days (adjust for production length)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -o md-%j.out
#SBATCH -e md-%j.err

module purge
module load gromacs/2024.3

echo "Start: $(date) | Node: $(hostname) | GPU: $CUDA_VISIBLE_DEVICES"

cd $SLURM_SUBMIT_DIR/GMX_PROLIG_MD
bash localrun.sh

echo "End: $(date)"
```

---

## PMF Workflow

PMF calculations have three stages. You can submit them as separate jobs with dependencies, or run them sequentially in one script.

### Option A: Single Script (Sequential)

```bash title="submit_pmf_all.sh"
#!/bin/bash
#SBATCH -J PRISM_PMF
#SBATCH -p <your_partition>
#SBATCH --time=168:00:00              # 7 days for all stages
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -o pmf-%j.out
#SBATCH -e pmf-%j.err

module purge
module load gromacs/2024.3

echo "Start: $(date) | Node: $(hostname) | GPU: $CUDA_VISIBLE_DEVICES"

cd $SLURM_SUBMIT_DIR/GMX_PROLIG_MD

# Stage 1: Steered MD
echo "=== Running Steered MD ==="
bash run_smd.sh

# Stage 2: Umbrella Sampling (sequential)
echo "=== Running Umbrella Sampling ==="
bash run_umbrella.sh

# Stage 3: WHAM Analysis
echo "=== Running WHAM Analysis ==="
bash run_wham.sh

echo "End: $(date)"
```

### Option B: Array Job for Umbrella Windows (Parallel)

For maximum throughput, submit each umbrella window as a separate GPU job:

```bash title="submit_umbrella_array.sh"
#!/bin/bash
#SBATCH -J PMF_umbrella
#SBATCH -p <your_partition>
#SBATCH --array=0-39                  # adjust to your number of windows
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH -o umbrella-%A_%a.out
#SBATCH -e umbrella-%A_%a.err

module purge
module load gromacs/2024.3

echo "Window ${SLURM_ARRAY_TASK_ID} | Node: $(hostname) | GPU: $CUDA_VISIBLE_DEVICES"

cd $SLURM_SUBMIT_DIR/GMX_PROLIG_MD/umbrella/window_${SLURM_ARRAY_TASK_ID}
gmx mdrun -deffnm umbrella -v -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} \
    -nb gpu -bonded gpu -pme gpu

echo "Window ${SLURM_ARRAY_TASK_ID} done: $(date)"
```

Workflow with job dependencies:

```bash
# Step 1: Run SMD first
JOB1=$(sbatch --parsable submit_smd.sh)

# Step 2: Umbrella array starts after SMD finishes
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 submit_umbrella_array.sh)

# Step 3: WHAM starts after all umbrella windows finish
sbatch --dependency=afterok:$JOB2 submit_wham.sh
```

---

## REST2 Replica Exchange

REST2 uses multiple replicas running simultaneously via `gmx mdrun -multidir`. This requires **one GPU per replica** or time-sharing on fewer GPUs.

```bash title="submit_rest2.sh"
#!/bin/bash
#SBATCH -J PRISM_REST2
#SBATCH -p <your_partition>
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16                   # one task per replica (match --replica-number)
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:4                  # adjust to available GPUs
#SBATCH --mem=64G
#SBATCH -o rest2-%j.out
#SBATCH -e rest2-%j.err

module purge
module load gromacs/2024.3

echo "Start: $(date) | Node: $(hostname)"
echo "Tasks: $SLURM_NTASKS | GPUs: $CUDA_VISIBLE_DEVICES"

WORKDIR=$SLURM_SUBMIT_DIR/GMX_PROLIG_MD
NREP=16                               # must match --replica-number

# Build multidir argument list
MULTIDIR=""
for i in $(seq 0 $((NREP - 1))); do
    MULTIDIR="$MULTIDIR $WORKDIR/replica_$i"
done

# Run REST2 with Hamiltonian replica exchange
gmx mdrun -deffnm md -multidir $MULTIDIR \
    -hrex -replex 1000 \
    -ntmpi 1 -ntomp ${SLURM_CPUS_PER_TASK} \
    -nb gpu -bonded gpu -pme gpu

echo "End: $(date)"
```

!!! note "GPU Allocation"
    REST2 with 16 replicas ideally needs 16 GPUs. If your cluster has fewer GPUs per node, either request multiple nodes or reduce the replica count. GROMACS will time-share replicas across available GPUs automatically.

---

## MM/PBSA

```bash title="submit_mmpbsa.sh"
#!/bin/bash
#SBATCH -J PRISM_MMPBSA
#SBATCH -p <your_partition>
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH -o mmpbsa-%j.out
#SBATCH -e mmpbsa-%j.err

module purge
module load gromacs/2024.3

echo "Start: $(date) | Node: $(hostname) | GPU: $CUDA_VISIBLE_DEVICES"

cd $SLURM_SUBMIT_DIR/GMX_PROLIG_MMPBSA
bash localrun.sh

echo "End: $(date)"
```

---

## Customization Guide

### Key Parameters to Adjust

| Parameter | What to change | How to decide |
| --- | --- | --- |
| `-p <partition>` | Your cluster's GPU partition name | Run `sinfo` to list available partitions |
| `--time` | Walltime limit | Standard MD: 3-7 days; PMF: 7-14 days |
| `--cpus-per-task` | CPU cores per task | 8-10 is typical for GPU jobs |
| `--gres=gpu:1` | GPU count and type | Use `gpu:a100:1` to request specific GPU type |
| `--mem` | Memory per node | 16G for standard MD; 64G+ for REST2 |
| `module load gromacs/...` | GROMACS version on your cluster | Run `module avail gromacs` to check |

### Monitoring Jobs

```bash
# Check job status
squeue -u $USER

# Check job output in real time
tail -f slurm-<jobid>.out

# Cancel a job
scancel <jobid>

# Check job efficiency after completion
seff <jobid>
```

### Email Notifications

Add these lines to receive email when jobs start, end, or fail:

```bash
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=your_email@example.com
```
