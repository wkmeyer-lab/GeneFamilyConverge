# Data/raw/

Immutable input files. **Never modified by the pipeline.**

## proteomes/

One FASTA per species. Naming: `<Species_name>.fasta`

Source: Ensembl BioMart (pilot). Document download URLs below as you add species.

| Species          | Source   | URL / accession              | Date downloaded |
|------------------|----------|------------------------------|-----------------|
| (add as needed)  |          |                              |                 |

## trees/

Rooted ultrametric species tree in Newick format. CAFE-5 requires ultrametric.

| File             | Source   | Notes                        |
|------------------|----------|------------------------------|
| (add as needed)  |          |                              |

## phenotypes/

Species-level phenotype tables (TSV). Columns must include species name
matching FASTA filenames and at least one phenotype column.

## metadata/

Any supplementary metadata files (literature tables, species lists from
VGP, etc.) that informed the phenotype assignments or species selection.

## Verification

After populating this directory, configure the run and preview the pipeline
without submitting anything:
```bash
convgeno init          # point it at your raw proteomes here
convgeno run -n        # dry run: prints the planned DAG (clean → OrthoFinder → CAFE-5)
```
