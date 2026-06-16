# Src/Reu — Reusable library code

Code here has **no hardcoded paths** and **no project-specific knowledge**.
Everything is passed as arguments. This code could be used by a different
project with different species, tools, or phenotypes.

## Python: `Src/Reu/python/`

Installable package (`convgeno`). Install in editable mode:

```bash
pip install -e Src/Reu/python/
```

Then import anywhere:

```python
from convgeno.io.fasta import parse_fasta
from convgeno.external.orthofinder import build_orthofinder_command
```

## R: `Src/Reu/r/`

Utility functions sourced by Loc R scripts:

```r
source("Src/Reu/r/read_cafe_outputs.R")
```

Kept simple intentionally — not a full R package. If the R code grows
substantially, it can be restructured into a package later.
