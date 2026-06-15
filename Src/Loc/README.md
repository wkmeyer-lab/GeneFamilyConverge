# Src/Loc — Project-specific pipeline code

This directory contains code specific to **this** pipeline run at Lehigh.
Scripts here are allowed to reference project paths, config files, and
the specific chaining of OrthoFinder → CAFE-5 / BadiRate → R analysis.

## Contents

- `configs/` — YAML config files for pipeline runs and tool paths
- `scripts/` — One script per pipeline step (Python upstream, R downstream)
- `workflows/` — Shell scripts and SLURM job files for local and HPC execution
- `reports/` — RMarkdown templates for pipeline result reports
