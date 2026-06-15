# Data/

All pipeline data lives here, separated into stages.

## Layout

| Directory     | Contents                                    | Git-tracked? |
|---------------|---------------------------------------------|--------------|
| `raw/`        | Immutable inputs. Never modified by pipeline | No (large)   |
| `interim/`    | Intermediate files between pipeline steps    | No           |
| `processed/`  | Final tool outputs and analysis results      | No (large)   |
| `logs/`       | Stdout/stderr from tool runs                 | No           |
| `manifests/`  | Checksums, tool versions, run metadata       | **Yes**      |

## Reproducing from scratch

1. Populate `raw/` using download instructions in `raw/README.md`
2. Verify checksums against `manifests/checksums.tsv`
3. Run the pipeline (see top-level README)
4. Compare new `manifests/` against committed versions
