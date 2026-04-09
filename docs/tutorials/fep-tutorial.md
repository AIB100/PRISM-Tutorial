# FEP Workflow Tutorial

**Level:** Advanced
**Time:** 2-4 hours (+ computation time)
**Topics:** Relative binding free energy, FEP, lambda windows, multi-estimator analysis

Learn how to calculate relative binding free energies between similar ligands using PRISM's complete FEP workflow.

## Objectives

In this tutorial, you will learn:

- How to prepare ligands for FEP calculations
- Setting up atom mapping and hybrid topology
- Running lambda window simulations
- Analyzing results with BAR/MBAR/TI estimators
- Validating convergence and quality

## Prerequisites

- Completed [Basic Tutorial](basic-tutorial.md) or familiar with PRISM basics
- GROMACS installed and configured
- 4+ GPUs available for parallel lambda window simulations
- Basic understanding of free energy calculations

## Background: The 42-38 System

We'll use the **42-38 ligand pair** from HIF-2α (Hypoxia-Inducible Factor 2α):

- **Ligand 42**: Reference ligand with a methoxyethylphenyl group
- **Ligand 38**: Mutant ligand with an ethylphenyl group (missing methoxy)

This is a **single-point mutation** ideal for FEP:
- Small structural change (removal of -OCH₃)
- High chemical similarity
- Well-defined binding mode

**Expected ΔG**: ~1-3 kcal/mol (mutant weaker binding)

## Step 1: Obtain Test System

Download the tutorial data:

```bash
# Clone PRISM repository (includes test data)
git clone https://github.com/AIB001/PRISM.git
cd PRISM/tests/gxf/FEP/unit_test/42-38

# Or download directly
wget https://github.com/AIB001/PRISM/raw/main/tests/gxf/FEP/unit_test/42-38/input/receptor.pdb
wget https://github.com/AIB001/PRISM/raw/main/tests/gxf/FEP/unit_test/42-38/input/42.pdb
wget https://github.com/AIB001/PRISM/raw/main/tests/gxf/FEP/unit_test/42-38/input/38.pdb
```

**Directory structure**:

```
42-38/
├── input/
│   ├── receptor.pdb    # HIF-2α protein (244 residues)
│   ├── 42.pdb          # Reference ligand (with -OCH₃)
│   └── 38.pdb          # Mutant ligand (without -OCH₃)
└── configs/
    ├── fep_gaff2.yaml  # FEP configuration
    └── case_gaff2.yaml # System configuration
```

## Step 2: Examine FEP Configuration

Let's look at the FEP configuration file:

```bash
cat configs/fep_gaff2.yaml
```

**Key parameters**:

```yaml
# Atom mapping
mapping:
  dist_cutoff: 0.6        # Distance threshold (Å)
  charge_cutoff: 0.05     # Charge difference threshold
  charge_common: mean     # Average charges for common atoms

# Lambda windows
lambda:
  strategy: decoupled     # Separate coulomb + vdw
  windows: 32             # Total windows (high accuracy)
  coul_windows: 12        # Electrostatic windows
  vdw_windows: 20         # Van der Waals windows

# Simulation
simulation:
  temperature: 310        # physiological temperature
  production_time_ns: 2.0 # Production per window
  dt: 0.002              # 2 fs timestep
```

**Why 32 windows?**
- More windows = better overlap = lower uncertainty
- 32 windows provide ~0.1 kcal/mol precision
- Trade-off: 3-4x more computation time vs 11 windows

## Step 3: Build FEP System

Create the FEP scaffold:

```bash
# From 42-38 directory
prism input/receptor.pdb input/42.pdb -o amber14sb_ol15-mut_gaff2 \
  --fep \
  --mutant input/38.pdb \
  --ligand-forcefield gaff2 \
  --forcefield amber14sb_OL15 \
  --fep-config configs/fep_gaff2.yaml \
  --config configs/case_gaff2.yaml
```

**What PRISM does**:

1. **Reads ligands**: Parses 42.pdb and 38.pdb
2. **Performs mapping**: Matches atoms by distance + charge
3. **Creates hybrid topology**: Single ITP with A/B states
4. **Sets up lambda windows**: Creates 32 window directories
5. **Generates scripts**: Automated run scripts for all legs

**Expected output** (takes 2-5 minutes):

```
✓ Atom mapping complete: 28 common, 4 transformed, 0 surrounding
✓ Hybrid topology generated
✓ FEP scaffold created: amber14sb_ol15-mut_gaff2/GMX_PROLIG_FEP/
  - bound/repeat1/ with 32 lambda windows
  - unbound/repeat1/ with 32 lambda windows
  - common/hybrid/ with topology files
```

## Step 4: Inspect Atom Mapping

Open the mapping visualization:

```bash
# View in browser
firefox amber14sb_ol15-mut_gaff2/GMX_PROLIG_FEP/common/hybrid/mapping.html
```

**What to check**:

1. **No gray atoms**: All atoms should be colored (no `classification: "unknown"`)
2. **Mapping statistics**:
   ```
   Ligand A: common=28, transformed=2, surrounding=0
   Ligand B: common=28, transformed=2, surrounding=0
   ```
3. **Chemical correctness**:
   - Common atoms: C, H, N, O backbone
   - Transformed A: O, CH₃ (methoxy group)
   - Transformed B: H (replacing methoxy)
4. **Total charge**: Should be ≈0 for neutral ligands

**Troubleshooting**:
- ❌ Gray atoms: Use MOL2 files instead of PDB (includes bond orders)
- ❌ Too many transformed atoms: Adjust `dist_cutoff` in config
- ❌ Wrong charges: Check ligand protonation states

## Step 5: Understand Directory Structure

Examine the created scaffold:

```bash
cd amber14sb_ol15-mut_gaff2/GMX_PROLIG_FEP

# Bound leg structure
ls bound/repeat1/
# window_00/  window_01/  ...  window_31/  (32 lambda windows)

# Unbound leg structure
ls unbound/repeat1/
# window_00/  window_01/  ...  window_31/

# Common hybrid files
ls common/hybrid/
# hybrid.itp  hybrid.gro  mapping.html
```

**Each lambda window contains**:

```
window_XX/
├── grompp.mdp        # Lambda-specific MDP parameters
├── topol.tpr         # GROMACS topology
├── conf.gro          # Initial coordinates
└── run.sh            # Execution script
```

**Lambda progression** (decoupled strategy):

- Windows 0-11: Coulomb transformation (λ: 0.0 → 1.0)
- Windows 12-31: VDW transformation (λ: 0.0 → 1.0)

## Step 6: Quick Test (100 Steps)

Before full production, test with short runs:

```bash
cd amber14sb_ol15-mut_gaff2/GMX_PROLIG_FEP

# Test bound leg, first window only
PRISM_MDRUN_NSTEPS=100 bash run_fep.sh bound 1

# Check if it completed
ls bound/repeat1/window_00/md.log
```

**Expected output**:

```
GROMACS reminder: ...
Running EM...
EM completed: 100 steps, max force = XXX kJ/mol/nm
Running NVT...
NVT completed: 100 steps
Running NPT...
NPT completed: 100 steps
Running production...
Production completed: 100 steps
```

**If this works**, your system is properly set up!

## Step 7: Run Full Production

For the complete calculation:

```bash
# Option 1: Run all (bound + unbound, all windows)
bash run_fep.sh all

# Option 2: Run bound leg only
bash run_fep.sh bound

# Option 3: Run specific repeat and window
bash run_fep.sh bound 1 5  # Bound, repeat1, window5
```

**Parallel execution**:

The script automatically runs 4 windows in parallel (configurable):

```bash
# Edit configs/fep_gaff2.yaml to change parallelization
execution:
  parallel_windows: 4   # Number of concurrent windows
  total_cpus: 56        # Total CPU cores
  num_gpus: 4           # Number of GPUs
```

**Estimated time** (32 windows, 2 ns each):

- Per window: 4-8 hours (depends on system size)
- Total: 32 windows × 6 hours ÷ 4 GPUs = **~48 hours**

## Step 8: Monitor Progress

Check simulation status:

```bash
# Check completion status
for i in {00..31}; do
  if [ -f bound/repeat1/window_${i}/production.log ]; then
    echo "Window $i: $(tail -1 bound/repeat1/window_${i}/production.log)"
  fi
done

# Or use the monitoring script
bash check_em_nvt_progress.sh amber14sb_ol15-mut_gaff2
```

**Look for**:
- ✅ "Finished mdrun" in production.log
- ✅ Reasonable EM convergence (max force < 1000 kJ/mol/nm)
- ❌ "Segmentation fault" or domain decomposition errors

## Step 9: Analyze Results

After simulations complete, run the analysis:

```bash
# Generate dhdl.xvg files first (if not auto-generated)
cd amber14sb_ol15-mut_gaff2/GMX_PROLIG_FEP

# Run multi-estimator analysis
python -m prism.fep.analysis.cli \
  --bound bound \
  --unbound unbound \
  --estimator BAR MBAR TI \
  --n-bootstrap 1000 \
  --n-jobs 8 \
  --output fep_results.html
```

**Analysis steps**:

1. **Parse dhdl.xvg**: Extract ∂H/∂λ from each window
2. **Calculate ΔG**: Apply BAR, MBAR, TI estimators
3. **Bootstrap**: Estimate uncertainties with 1000 resamples
4. **Generate report**: Create interactive HTML

**Expected results**:

```
=== FEP Analysis Results ===

Bound leg ΔG:   -5.23 ± 0.15 kcal/mol
Unbound leg ΔG: -3.45 ± 0.12 kcal/mol

ΔΔG (binding): 1.78 ± 0.20 kcal/mol

Estimator breakdown:
  BAR:  1.75 ± 0.18 kcal/mol
  MBAR: 1.78 ± 0.20 kcal/mol
  TI:   1.82 ± 0.22 kcal/mol
```

**Interpretation**:

- Positive ΔΔG: Ligand 38 binds weaker than ligand 42
- Magnitude: ~1.8 kcal/mol (small but measurable)
- Consistency: All 3 estimators agree (good sign!)

## Step 10: Validate Quality

Open the HTML report:

```bash
firefox fep_results.html
```

**Quality checks**:

1. **Overlap matrix**: Green/yellow cells = good overlap
2. **dhdl plots**: Smooth curves without jumps
3. **Convergence**: ΔG stabilizes over time
4. **Uncertainties**: SE < 1.0 kcal/mol (acceptable)

**Warning signs**:

- ⚠️ Red cells in overlap matrix → Add more lambda windows
- ⚠️ Large error bars (> 1 kcal/mol) → Increase production time
- ⚠️ Diverging ΔG over time → Check simulation stability

## Step 11: Compare with Experiment

If experimental data is available:

```
Experimental ΔΔG (42→38): 2.1 ± 0.3 kcal/mol
Calculated ΔΔG (42→38):  1.8 ± 0.2 kcal/mol

Difference: 0.3 kcal/mol (within error bars)
```

**Good agreement**! The FEP calculation correctly predicts:
- Ligand 38 binds weaker
- Magnitude is reasonable
- Uncertainty captures experimental error

## Advanced: Force Field Comparison

Test how force field choice affects results:

```bash
# GAFF2 result (from above)
ΔΔG_GAFF2 = 1.8 ± 0.2 kcal/mol

# Now try OpenFF
prism input/receptor.pdb input/42.pdb -o amber14sb_ol15-mut_openff \
  --fep --mutant input/38.pdb \
  --ligand-forcefield openff \
  --forcefield amber14sb_OL15 \
  --fep-config configs/fep_openff.yaml

# OpenFF result
ΔΔG_OpenFF = 2.1 ± 0.3 kcal/mol

# Compare
GAFF2: 1.8 ± 0.2 kcal/mol
OpenFF: 2.1 ± 0.3 kcal/mol
Experiment: 2.1 ± 0.3 kcal/mol

→ OpenFF closer to experiment for this system!
```

## Troubleshooting

### Issue: Mapping shows gray atoms

**Symptom**: `classification: "unknown"` in mapping HTML

**Cause**: Missing bond order information (PDB files)

**Fix**:
```bash
# Use MOL2 files instead
obabel input/42.pdb -O input/42.mol2 -h
obabel input/38.pdb -O input/38.mol2 -h

# Rebuild with MOL2 files
prism input/receptor.pdb input/42.mol2 -o output \
  --fep --mutant input/38.mol2 \
  --ligand-forcefield gaff2
```

### Issue: Simulation crashes during EM

**Symptom**: "Segmentation fault" or nan coordinates

**Cause**: Incorrect hybrid topology or bad starting structure

**Fix**:
```bash
# Check hybrid topology
grep "atoms" common/hybrid/hybrid.itp -A 50

# Verify B-state masses are present
# Every transformed atom should have mass_b > 0

# If missing masses, rebuild with latest PRISM
```

### Issue: Poor overlap between windows

**Symptom**: Red cells in overlap matrix

**Cause**: Too few lambda windows or poor lambda spacing

**Fix**:
```yaml
# Add more windows in fep.yaml
lambda:
  windows: 21        # Increase from 11
  coul_windows: 8    # More coulomb windows
  vdw_windows: 13     # More vdw windows
```

### Issue: High uncertainty in ΔG

**Symptom**: Standard error > 1.0 kcal/mol

**Cause**: Insufficient sampling per window

**Fix**:
```yaml
# Increase production time
simulation:
  production_time_ns: 5.0  # Increase from 2.0
```

## Summary

In this tutorial, you:

✅ Built a complete FEP system for the 42-38 ligand pair
✅ Inspected atom mapping and hybrid topology
✅ Ran lambda window simulations (or learned how to)
✅ Analyzed results with BAR/MBAR/TI estimators
✅ Validated convergence and quality
✅ Compared calculated ΔΔG with experiment

**Key takeaways**:

- FEP is **powerful** for relative binding free energies
- **Proper setup** is critical (mapping, hybrid topology)
- **Quality checks** ensure reliable results
- **Multiple estimators** provide robustness
- **Force field choice** affects accuracy

## Next Steps

- Try the [Force Field Tutorial](force-field-tutorial.md) to compare multiple force fields
- Learn about [PMF Calculations](pmf-tutorial.md) for absolute binding free energy
- Explore [Batch Processing](batch-tutorial.md) for high-throughput FEP
- Study the [API Reference](../api/fep.md) for custom FEP workflows

## Additional Resources

- [FEP User Guide](../user-guide/fep-calculations.md) - Complete FEP documentation
- [Test Systems](https://github.com/AIB001/PRISM/tree/main/tests/gxf/FEP/unit_test) - More examples
- [Analysis Tools](../user-guide/analysis-tools.md) - Post-processing FEP data

<div class="whats-next" markdown>

## What's Next

- [Learn about PMF Calculations for absolute binding free energy](pmf-tutorial.md)
- [Compare force fields in the Force Field Tutorial](force-field-tutorial.md)
- [Explore advanced analysis techniques](../user-guide/analysis-tools.md)

</div>
