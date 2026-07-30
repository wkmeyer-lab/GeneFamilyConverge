# Species tree

The pipeline analyzes gene-family change on an **ultrametric** species tree: a
rooted tree in which every tip is the same total distance from the root, with
branch lengths in units of time.

You have three options for supplying this tree, and you pick one while setting up
the run with `convgeno init`.

| Option | What you provide | Branch-length units | r8s |
|--------|------------------|---------------------|-----|
| 1 | Nothing (the default) | — | Yes |
| 2 | A non-ultrametric tree | Substitutions per site | Yes |
| 3 | An ultrametric tree | Millions of years | No |

Your tree must be in **Newick** format, and its tip labels must exactly match your
proteome file names (without the file extension).

`convgeno init` sets the chosen option up for you from your answers. Each option
below also lists the flags the species-tree step uses; you can pass these directly
if you run that step yourself.

## Option 1 — use the OrthoFinder species tree

OrthoFinder builds the species tree from your proteomes, and the pipeline runs r8s
to make it ultrametric.

- **At `convgeno init`:** press Enter at the species-tree prompt.
- **Flags:** none (this is the default).

## Option 2 — provide your own non-ultrametric tree

Provide a species tree whose branch lengths are in **substitutions per site** (the
branch-length units produced by OrthoFinder and most tree-building programs). The
pipeline runs r8s to convert it into an ultrametric time tree.

You may also provide the **number of sites**: the number of columns (sites) in the
alignment that your tree's branch lengths were estimated from. If you do not provide
it, the pipeline uses the number of sites from the OrthoFinder alignment and runs r8s.

- **At `convgeno init`:** give the path to your tree, answer **no** when asked
  whether the tree is ultrametric, then enter the number of sites (or press Enter to
  leave it unset).
- **Flags:** `--input-tree <path>`, and optionally `--nsites <N>`.

## Option 3 — provide your own ultrametric tree

Provide a species tree whose branch lengths are in **millions of years** (an
ultrametric tree). The pipeline uses it directly and does not run r8s.

- **At `convgeno init`:** give the path to your tree and answer **yes** when asked
  whether the tree is ultrametric.
- **Flags:** `--input-tree <path> --assume-ultrametric`.

If you mark a tree as ultrametric but it is not, the pipeline reports this and dates
it with r8s instead (as in Option 2).
