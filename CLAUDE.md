# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Pipeline for characterizing gene family copy number across species and associating
copy-number variation with convergent phenotypes (Aim 1a of the Meyer Lab NSF CAREER
project). The scientific flow is: proteomes → filter isoforms → OrthoFinder (gene
families) → CAFE-5 / BadiRate (turnover modeling) → R analysis/figures. The downstream
R association step (all stubs today) is intended to link copy-number variation to
convergent phenotypes; the planned association method is RERconverge-adapted branch
association.

**Implementation status matters when reading this repo.** Only part of that flow is
built. The `convgeno` Python package is the implemented, heavily-tested core, and two
things in it actually work end-to-end:

1. **Longest-isoform filtering** — `convgeno.io.fasta` (~720 lines, large test suite).
   This is the only *scientific* pipeline step that is complete. It does single-pass
   and two-pass low-memory (`memory_mode='low'`, `_filter_low_memory`) filtering,
   transparent gzip/bz2 decompression via magic-byte sniffing, Ensembl/NCBI header
   auto-detection, and multi-file-per-species batch processing.
2. **OrthoFinder SLURM orchestration** — generating and submitting single- and
   multi-node SLURM jobs (the `slurm/` package + `orthofinder` CLI). This is plumbing,
   not a scientific step.

Almost everything else is a **stub** — docstring-only "functions to implement" files or
scripts that print "not yet implemented" and exit 1. This includes the
`convgeno.external.cafe` / `convgeno.external.badirate` modules, CAFE/BadiRate/most
`Src/Loc/scripts/*.py`, all of `Src/Reu/r/`, and several `convgeno` support modules
(see the module notes below for which are real vs. stub). Don't assume a step works
because the README diagram or `run_local.sh` references it; check the file first.

## Commands

Environment setup:
```bash
conda env create -f environment.yml
conda activate convgeno
pip install -e Src/Reu/python/          # editable install; also done by environment.yml
```

Python package (all commands run from `Src/Reu/python/`):
```bash
pytest                                   # full suite (361 tests), config in pyproject.toml
pytest tests/test_slurm_discovery.py     # single test file
pytest tests/test_runtime.py::test_name  # single test
pytest --cov=convgeno                     # with coverage
ruff check .                              # lint (rules E,F,W,I,N,UP,B,SIM; line-length 88)
ruff format .                             # format
```

The `convgeno` CLI (installed as a console script) — `init` + OrthoFinder SLURM only:
```bash
convgeno init                            # interactive wizard → writes pipeline_config.yaml
convgeno orthofinder generate            # write SLURM script(s), do not submit
convgeno orthofinder run                 # write + sbatch submit (prompts unless -y)
convgeno orthofinder generate --multinode   # 3-script prepare→search-array→resume chain
```

Isoform filtering has **no `convgeno` subcommand**; its CLI is a Loc script that imports
`convgeno.io.fasta` (run from the repo root, after `pip install -e Src/Reu/python/`):
```bash
# batch: one FASTA per species in a directory, optional per-species stats JSON
python Src/Loc/scripts/filter_isoforms.py dir <input_dir> <output_dir> [--stats-json stats.json]
# single species: one or more input files (e.g. per-chromosome) → one output
python Src/Loc/scripts/filter_isoforms.py files -i in1.fa[.gz] [in2.fa …] -o out.fa
# shared flags: --format {auto,ensembl,ncbi}  --on-duplicate {error,warn,skip}  --log-level  --log-file
```

## Architecture

### Reu vs Loc — the core convention
- `Src/Reu/` = **reusable** library code with **no hardcoded paths** and no
  project-specific knowledge; everything is passed as arguments. Contains the
  installable `convgeno` Python package (`Src/Reu/python/`) and standalone R utility
  scripts (`Src/Reu/r/`, sourced, not a package).
- `Src/Loc/` = **project-specific** glue: this run's configs, per-step scripts,
  SLURM/shell workflows, and RMarkdown reports. Allowed to reference concrete paths and
  chain tools together.
- Python is upstream (ingestion, filtering, format prep, tool invocation); R is
  downstream (parsing results, stats, figures). **The two languages communicate only
  via files** (TSV, FASTA, Newick) — never in-process.

### `convgeno` package layout (`Src/Reu/python/src/convgeno/`)
- `cli/` — argparse entry point (`cli/__init__.py:main`) dispatching to `init_cmd.py`
  (interactive wizard) and `orthofinder_cmd.py` (generate/submit handlers).
- `slurm/` — the architectural heart:
  - `config.py` — frozen dataclasses `SlurmConfig` / `PipelineConfig` with YAML
    load/save. `mem` and `mem_per_cpu` are mutually exclusive (enforced in `__post_init__`).
  - `discovery.py` — probes the cluster via `sinfo` **and `scontrol`**: `sinfo` for
    partitions, per-node CPUs, and memory limits; `scontrol show node` to read
    `ThreadsPerCore` for the physical-vs-logical core derivation on hyperthreaded nodes.
  - `runtime.py` — **single source of truth for conda-activation shell code.** Compute
    nodes start bare (no conda on PATH), so `$(conda info --base)` is a bootstrap
    paradox. Instead absolute conda paths are detected at `init` time and embedded
    literally into every generated script. No other module should contain inline conda
    bootstrap logic.
  - `script_generator.py` — single-node OrthoFinder SLURM script.
  - `multinode_generator.py` — the 3-job chain (prepare `-op` → search array →
    resume `-b`). Contains carefully-tuned regexes that filter real DIAMOND/BLAST search
    commands out of OrthoFinder's log prose; changing them can silently break the array.
- `external/` — tool config + wrappers. `config.py` (`OrthoFinderConfig`) is real;
  `cafe.py`, `badirate.py`, `orthofinder.py` are planned stubs. Note `orthofinder.py` is
  a stub *by design* — actual execution goes through the generated SLURM scripts, not an
  in-process wrapper.
- Support modules — **verify before trusting**, several are stubs like `Src/Loc/scripts`:
  - **Real:** `io/fasta.py` (isoform filtering, above), `utils/logging.py` (`setup_logging`),
    `validation/orthofinder_inputs.py` (pre-submission input checks used by
    `orthofinder_cmd.py`).
  - **Stub (docstring-only):** `io/tables.py`, `io/paths.py`, `validation/proteomes.py`,
    `validation/trees.py`, `manifests/checksums.py`, `manifests/run_metadata.py`,
    `utils/command_runner.py`, `utils/config.py`.

### SLURM job submission behavior (`orthofinder_cmd.py`)
- `generate` writes scripts; `run` also submits via `sbatch`. Default flow prints the
  script and asks for confirmation; `-y/--yes` skips it.
- **Automatic account fallback:** if a script has `#SBATCH --account=…` and SLURM
  rejects it, submission retries once with the account line stripped. In multi-node mode,
  once the first job succeeds without an account, the rest of the dependency chain is
  submitted without it too, to stay consistent.
- OrthoFinder refuses a pre-existing non-default `-o` output dir, so `init` defaults to a
  fresh timestamped `Data/processed/orthofinder_single_<timestamp>/` and
  `validate_orthofinder_output_dir` rejects existing dirs and paths with quote chars.

### Data/ stages
`raw/` (immutable inputs) → `interim/` → `processed/` → plus `logs/` and `manifests/`.
Only `manifests/` (checksums, tool versions, run metadata) is git-tracked; raw/interim/
processed are gitignored (large/generated), preserved by `.gitkeep` placeholders.

## Conventions
- Config dataclasses ignore unknown YAML keys in `from_dict` for forward compatibility;
  `to_dict`/`save` round-trip through YAML.
- Optional SLURM account is normalized (`normalize_optional_account`): blank, `"null"`,
  `"none"`, `None` all mean "emit no `--account` line."
- Tests use pytest `tmp_path` fixtures (see `tests/conftest.py` for sample Ensembl/NCBI
  FASTA fixtures); no test writes outside its tmp dir.
- Pipeline-step scripts use **descriptive names, not numbered prefixes** (e.g.
  `filter_isoforms.py`, not `step01_filter_isoforms.py`).
- `*.local.yaml` and `tool_paths.local.yaml` are gitignored for local overrides.
