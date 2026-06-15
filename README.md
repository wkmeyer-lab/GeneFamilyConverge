# Gene Family Copy Number Pipeline

Pipeline for characterizing gene family copy number across species and
associating variation in copy number with convergent phenotypes.

Supports **Aim 1a** of the Meyer Lab NSF CAREER project.

## Pipeline overview

```
Proteomes → filter isoforms → OrthoFinder → gene family counts
                                                  ↓
                                           ┌──────┴──────┐
                                           ↓              ↓
                                        CAFE-5        BadiRate
                                     (size/turnover)  (model comparison)
                                           ↓              ↓
                                           └──────┬──────┘
                                                  ↓
                                        R analysis + figures
```

## Confirmed tools

| Tool        | Purpose                                         |
|-------------|-------------------------------------------------|
| OrthoFinder | Define gene families via ortholog clustering     |
| CAFE-5      | Model gene family size/turnover on a phylogeny   |
| BadiRate    | Compare null vs phenotype-driven evolution models |

## Language split

- **Python**: upstream (ingestion, filtering, format prep, CLI tool invocation)
- **R**: downstream (parse results, statistical analysis, visualization)
- Communication between languages is always via files (TSV, FASTA, Newick)

## Structure

```
Src/Loc/     Project-specific scripts, configs, workflows, reports
Src/Reu/     Reusable library code (Python package + R utility scripts)
Data/        raw → interim → processed, plus logs and manifests
```

## Setup

```bash
conda env create -f environment.yml
conda activate genefam-pipeline

# Install the reusable Python package in editable mode
pip install -e Src/Reu/python/

# R dependencies restored via renv on first R session
```
