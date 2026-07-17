# Calibrating the DIAMOND throughput constant

The multi-node search phase sizes each array task's `--time` from a per-core
DIAMOND rate, `throughput_const` (bytes² per core-second):

```
seconds = max_bucket_cost / (cores · throughput_const) · margin
```

`convgeno` ships a default calibrated on Lehigh Sol (`hawkcpu`), but DIAMOND's
speed varies by CPU, DIAMOND version, and sensitivity settings. **Recalibrate on
your own cluster** for accurate walltime requests — this kit measures the
constant empirically, and you then set it in your config.

## What it measures

For each sampled `diamond blastp` command the prepare phase emitted, the cost is
`bytes(Speciesᵢ.fa) · bytes(Speciesⱼ.fa)` (the exact model the generator uses).
The kit runs a random sample of commands **one at a time, single-threaded
(`-p 1`)**, times each, and reports:

```
throughput_const = Σ cost / Σ elapsed_seconds
```

plus the per-command min/median/max rates as a quality check.

## Prerequisites

1. Run the **prepare phase** so the emitted commands and the WorkingDirectory
   exist:
   ```bash
   convgeno orthofinder generate            # writes slurm_scripts/orthofinder_prepare.sh
   sbatch slurm_scripts/orthofinder_prepare.sh   # wait for COMPLETED
   ```
   This produces `<output_parent>/<run>_search_commands.txt` and the
   `WorkingDirectory` (with `Species*.fa`), referenced by
   `<output_parent>/<run>_working_dir_path.txt`.
2. The `convgeno` conda env (with `diamond`) must be activatable.

## Run it

### Option A — SLURM batch (recommended; 100 sequential searches take hours)

Edit `throughput_calibration.sbatch` (the `# EDIT:` lines: partition, and the
conda bootstrap — easiest is to paste the bootstrap block from the top of your
generated `slurm_scripts/orthofinder_prepare.sh`), then from the repo root:

```bash
sbatch tools/throughput-constant/throughput_calibration.sbatch
```

The wrapper resolves the prepare outputs from `pipeline_config.yaml` itself.

### Option B — interactive

```bash
salloc --partition=<your_partition> --cpus-per-task=2 --mem=8G --time=12:00:00
# then, inside the allocation, activate the convgeno env and:
OUT=$(python -c "import yaml;print(yaml.safe_load(open('pipeline_config.yaml'))['orthofinder']['output_dir'])")
PARENT=$(dirname "$OUT"); RUN=$(basename "$OUT")
bash tools/throughput-constant/throughput_calibration.sh \
    "$PARENT/${RUN}_search_commands.txt" \
    "$(cat "$PARENT/${RUN}_working_dir_path.txt")" \
    100
```

Run it on the **same partition** your search array will use — the number is
CPU-specific. `N=100` is a good sample; a smaller `N` is faster but noisier.

## Read the result

Look for:

```
>>> throughput_const = 6.055387e+11 bytes^2 per core-second
per-command rate min/median/max: 4.89e+11 / 6.08e+11 / 7.70e+11
```

**Sanity checks that the value is trustworthy** (see also the reliability notes
in `docs/orthofinder/search-phase.md`):

- **All sampled commands succeeded** — the estimate isn't biased by a
  failure-skewed subset.
- **Aggregate ≈ median** — if the cost-weighted aggregate and the per-command
  median are within a few percent, the estimate is representative rather than
  dominated by a few unusual commands.
- **A moderate min↔max spread is expected** — real runtime depends on more than
  FASTA byte size (sequence composition, hit counts, cache/filesystem effects).

## Apply it

Set the measured **aggregate** value in `pipeline_config.yaml`:

```yaml
multinode:
  throughput_const: 6.055387e11
```

The generator uses this instead of the built-in default (see
`derive_search_walltime`). No code edit is needed.

## Interpretation

This measures **isolated, single-core** throughput. When many single-threaded
DIAMOND searches run concurrently on a full node, the effective per-core rate can
drop (memory bandwidth, cache, filesystem contention). The walltime model keeps a
safety `margin` (default 1.5) for exactly this reason — do **not** shrink the
constant further to compensate unless real search-array runs show the margin is
insufficient. An over-estimate only costs queue time; an under-estimate risks a
task timeout (recoverable — the resume completeness gate reports the missing
searches and a resubmit line).
