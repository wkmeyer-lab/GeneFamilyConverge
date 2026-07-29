"""
convgeno.external.r8s
~~~~~~~~~~~~~~~~~~~~~~

Make an OrthoFinder species tree ultrametric with r8s, for input to CAFE-5.

CAFE-5 requires a rooted, binary, **ultrametric** time tree, but OrthoFinder
emits an *additive* (non-ultrametric) ``SpeciesTree_rooted.txt`` whose branch
lengths are expected substitutions per site.  r8s rescales those branch lengths
into time units via penalized likelihood, using one or more fossil/literature
calibration points.  This module reproduces (and hardens) the workflow from the
CAFE-5 tutorial's ``prep_r8s.py`` + ``r8s`` steps.

Capabilities
------------
- Derive ``nsites`` directly from the concatenated MSA that produced the tree
  (``count_alignment_sites``) — the number of *alignment columns*, which is what
  r8s's ``blformat lengths=persite`` needs, not a sum of raw sequence lengths.
- Locate OrthoFinder's species tree + species-tree alignment inside a
  ``Results_*`` directory (``find_results_dir`` / ``find_species_tree_files``).
- Strip FastTree support values / internal labels that r8s cannot parse
  (``sanitize_tree_for_r8s``).
- Generate the NEXUS-style r8s control file with ``fixage`` (point) and/or
  ``constrain`` (windowed) calibrations (``build_r8s_control``).
- Robustly extract the dated Newick from r8s output
  (``parse_ultrametric_tree``) — a resilient replacement for the tutorial's
  fragile ``tail -n 1 | cut -c 16-``.
- Orchestrate the whole step (``make_ultrametric``).

Dependencies
------------
- BioPython (``biopython >= 1.80``): ``Bio.Phylo`` for Newick parse/emit and
  ``Bio.SeqIO`` (via :func:`convgeno.io.fasta.open_fasta`) for alignment I/O.
- The r8s binary itself is **not** on conda; it must be installed from source /
  a prebuilt Linux binary and its path supplied by the caller.  Its invocation
  goes through :mod:`convgeno.utils.command_runner` (imported lazily so this
  module stays importable before that module is implemented).
"""

from __future__ import annotations

import logging
import re
from collections.abc import Iterable
from io import StringIO
from pathlib import Path
from typing import NamedTuple

from Bio import Phylo, SeqIO

from convgeno.io.fasta import open_fasta

logger = logging.getLogger(__name__)

__all__ = [
    "Calibration",
    "count_alignment_sites",
    "find_results_dir",
    "find_species_tree_files",
    "sanitize_tree_for_r8s",
    "root_spanning_taxa",
    "build_r8s_control",
    "build_command",
    "parse_ultrametric_tree",
    "make_ultrametric",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

#: Location of the rooted species tree relative to an OrthoFinder ``Results_*``.
SPECIES_TREE_RELPATH = Path("Species_Tree") / "SpeciesTree_rooted.txt"

#: Location of the concatenated species-tree MSA (MSA mode only) relative to a
#: ``Results_*`` directory.  Its column count *is* ``nsites``.
ALIGNMENT_RELPATH = Path("MultipleSequenceAlignments") / "SpeciesTreeAlignment.fa"

#: Name given to the tree inside the generated control file (matches prep_r8s.py).
_TREE_NAME = "nj_tree"

#: Matches an r8s ``tree <name> = <newick>;`` line and captures the Newick.
_TREE_LINE_RE = re.compile(r"tree\s+\S+\s*=\s*(.+;)", re.IGNORECASE)

#: Matches a leading rooting comment (``[&R]`` / ``[&U]``) on a Newick string.
_ROOTING_COMMENT_RE = re.compile(r"^\[&[RU]\]\s*")

# ---------------------------------------------------------------------------
# Calibration record
# ---------------------------------------------------------------------------


class Calibration(NamedTuple):
    """A single node calibration for r8s.

    The calibrated node is the most recent common ancestor (MRCA) of ``taxa``
    (two tip labels that must both exist in the tree).  Provide either a fixed
    ``age`` (→ r8s ``fixage``) or an age window via ``min_age`` / ``max_age``
    (→ r8s ``constrain``).  All ages are in the same time unit (e.g. Myr).
    """

    name: str
    taxa: tuple[str, str]
    age: float | None = None
    min_age: float | None = None
    max_age: float | None = None


# ===================================================================
#  1.  nsites — the number of ALIGNMENT COLUMNS
# ===================================================================


def count_alignment_sites(alignment_path: Path | str) -> int:
    """Return ``nsites``: the number of columns in a concatenated alignment.

    This is the value r8s needs for ``blformat lengths=persite nsites=N``.  It
    is the aligned **width** (number of columns) of the supermatrix that
    produced the tree's branch lengths — *not* a sum of unaligned residue
    counts, and *not* a per-file maximum.  In an alignment every record is the
    same width (gap-padded), so we read the first record's length and assert the
    rest match.

    Parameters
    ----------
    alignment_path : Path or str
        Path to ``SpeciesTreeAlignment.fa`` (OrthoFinder MSA mode).  Transparent
        gzip/bz2 decompression is handled by :func:`convgeno.io.fasta.open_fasta`.

    Returns
    -------
    int
        The number of alignment columns.

    Raises
    ------
    FileNotFoundError
        If the alignment file does not exist.
    ValueError
        If the file is empty or not rectangular (i.e. not a real alignment).
    """
    alignment_path = Path(alignment_path)
    if not alignment_path.is_file():
        raise FileNotFoundError(f"Species-tree alignment not found: {alignment_path}")

    logger.info("Counting alignment columns (nsites) in %s", alignment_path)

    width: int | None = None
    n_records = 0
    with open_fasta(alignment_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_len = len(record.seq)
            n_records += 1
            if width is None:
                width = seq_len
            elif seq_len != width:
                raise ValueError(
                    f"Alignment is not rectangular: record '{record.id}' has "
                    f"length {seq_len}, expected {width}. "
                    f"{alignment_path} does not look like an aligned FASTA."
                )

    if not n_records or width is None:
        raise ValueError(f"No sequences found in alignment: {alignment_path}")

    logger.info("nsites = %d (from %d aligned sequences)", width, n_records)
    return width


# ===================================================================
#  2.  Locating OrthoFinder output files
# ===================================================================


def _read_pointer_working_dirs(root: Path) -> list[Path]:
    """Read multi-node ``WorkingDirectory`` pointer(s) for *root*, if present.

    The multi-node prepare job writes ``<parent>/<name>_working_dir_path.txt``
    (a sibling of the output dir) holding the absolute path to OrthoFinder's
    ``WorkingDirectory``; the resume job runs ``-b`` in place there, so the final
    results live at ``<WorkingDirectory>/OrthoFinder/Results_*``.
    """
    working_dirs: list[Path] = []
    pointers = [
        root.parent / f"{root.name}_working_dir_path.txt",
        *root.glob("*_working_dir_path.txt"),
    ]
    for pointer in pointers:
        try:
            if pointer.is_file():
                wd = Path(pointer.read_text(encoding="utf-8").strip())
                if wd.is_dir():
                    working_dirs.append(wd)
        except OSError:
            continue
    return working_dirs


def _candidate_results_dirs(root: Path) -> list[Path]:
    """Return plausible ``Results_*`` directories for *root* (deduplicated).

    Covers both workflow layouts without a full recursive walk:

    - **single-node**: results at ``<root>/Results_*`` (OrthoFinder ``-o root``);
    - **multi-node**: ``-b`` runs in place, so results are at
      ``<WorkingDirectory>/OrthoFinder/Results_*`` (much deeper) — the
      WorkingDirectory is read from the sibling ``*_working_dir_path.txt``
      pointer.

    Also accepts being pointed straight at a ``Results_*`` directory.
    """
    candidates = [
        root,
        *sorted(root.glob("Results_*")),
        *sorted(root.glob("OrthoFinder/Results_*")),
    ]
    for wd in _read_pointer_working_dirs(root):
        candidates += [
            wd,
            *sorted(wd.glob("Results_*")),
            *sorted(wd.glob("OrthoFinder/Results_*")),
        ]
    seen: set[Path] = set()
    uniq: list[Path] = []
    for candidate in candidates:
        if candidate.is_dir() and candidate not in seen:
            seen.add(candidate)
            uniq.append(candidate)
    return uniq


def _recursive_results_dirs(root: Path) -> list[Path]:
    """Last-resort: any ``Results_*`` under *root* that holds a species tree.

    Handles arbitrarily deep layouts (e.g. the multi-node
    ``Results_*/WorkingDirectory/OrthoFinder/Results_*`` nesting) when the
    shallow + pointer checks miss. Bounded to *root*; only used as a fallback.
    """
    pattern = f"**/{SPECIES_TREE_RELPATH.as_posix()}"
    return [tree.parent.parent for tree in root.glob(pattern) if tree.is_file()]


def find_results_dir(orthofinder_output_dir: Path | str) -> Path:
    """Find the OrthoFinder ``Results_*`` directory holding the species tree.

    Robust to both the single-node (``<dir>/Results_*``) and multi-node
    (``<WorkingDirectory>/OrthoFinder/Results_*``, via the
    ``*_working_dir_path.txt`` pointer) layouts, and to being handed a
    ``Results_*`` directory directly. If several candidates contain a species
    tree, the most recently modified one wins. Fulfils the behaviour stubbed in
    :mod:`convgeno.external.orthofinder`.

    Raises
    ------
    NotADirectoryError
        If *orthofinder_output_dir* is not a directory.
    FileNotFoundError
        If no directory containing ``Species_Tree/SpeciesTree_rooted.txt`` is found.
    """
    root = Path(orthofinder_output_dir)
    if not root.is_dir():
        raise NotADirectoryError(f"OrthoFinder output directory not found: {root}")

    valid = [
        d for d in _candidate_results_dirs(root) if (d / SPECIES_TREE_RELPATH).is_file()
    ]
    if not valid:
        # Deep / unknown layout — walk (bounded to root) as a last resort.
        valid = _recursive_results_dirs(root)
    if not valid:
        raise FileNotFoundError(
            f"No OrthoFinder Results directory containing "
            f"{SPECIES_TREE_RELPATH.as_posix()} found under {root}. Checked the "
            "single-node layout (<dir>/Results_*), the multi-node WorkingDirectory "
            "pointer (*_working_dir_path.txt), and a recursive search. Has "
            "OrthoFinder finished and produced a species tree?"
        )

    results_dir = max(valid, key=lambda d: d.stat().st_mtime)
    if len(valid) > 1:
        logger.info(
            "Multiple OrthoFinder Results dirs found; using newest: %s", results_dir
        )
    return results_dir


def find_species_tree_files(orthofinder_output_dir: Path | str) -> tuple[Path, Path]:
    """Return ``(species_tree, species_tree_alignment)`` for an OrthoFinder run.

    The alignment is required so ``nsites`` can be derived; its absence almost
    always means OrthoFinder was run in the legacy distance-based mode
    (``-M dendroblast``) rather than MSA mode.

    Raises
    ------
    FileNotFoundError
        If no ``Results_*`` dir is found, or the run has no concatenated
        species-tree alignment.
    """
    results_dir = find_results_dir(orthofinder_output_dir)
    species_tree = results_dir / SPECIES_TREE_RELPATH
    alignment = results_dir / ALIGNMENT_RELPATH

    if not alignment.is_file():
        raise FileNotFoundError(
            f"Concatenated species-tree alignment not found: {alignment}\n"
            "This usually means OrthoFinder was run WITHOUT MSA mode (legacy "
            "'-M dendroblast'), so 'nsites' cannot be derived from an alignment. "
            "Re-run OrthoFinder with '-M msa', pass an explicit nsites, or use a "
            "distance-based ultrametric method instead of r8s."
        )

    logger.info("Species tree: %s", species_tree)
    logger.info("Species-tree alignment: %s", alignment)
    return species_tree, alignment


# ===================================================================
#  3.  Tree sanitisation
# ===================================================================


def _read_newick_text(source: Path | str) -> str:
    """Read Newick text from a path, or return *source* if it is Newick itself."""
    if isinstance(source, Path):
        return source.read_text(encoding="utf-8")
    text = str(source)
    if "(" in text or ";" in text:
        return text
    candidate = Path(text)
    if candidate.is_file():
        return candidate.read_text(encoding="utf-8")
    return text


def sanitize_tree_for_r8s(newick_in: Path | str) -> str:
    """Return a bare rooted Newick safe to feed to r8s.

    FastTree-based OrthoFinder species trees carry internal-node support values
    (parsed by ``Bio.Phylo`` as node names/confidences).  r8s cannot parse those
    and will error, so we strip every internal-node name/confidence/comment
    while preserving topology, tip labels, and all branch lengths.

    Parameters
    ----------
    newick_in : Path or str
        A path to a Newick file, or the Newick string itself.

    Returns
    -------
    str
        A single-line Newick ending in ``;`` (no ``[&R]`` prefix — that is added
        by :func:`build_r8s_control`).
    """
    text = _read_newick_text(newick_in)
    tree = Phylo.read(StringIO(text), "newick")

    for clade in tree.get_nonterminals():
        clade.name = None
        clade.confidence = None
        clade.comment = None
    for clade in tree.get_terminals():
        clade.confidence = None
        clade.comment = None

    out = StringIO()
    Phylo.write(tree, out, "newick", format_branch_length="%1.8f")
    newick = out.getvalue().strip()
    if not newick.endswith(";"):
        newick += ";"
    return newick


def _root_spanning_taxa(tree) -> tuple[str, str]:
    """Two tip labels whose MRCA is the tree root (one per side of the root)."""
    children = list(tree.root.clades)
    if len(children) < 2:
        raise ValueError(
            "Cannot anchor the root: the tree root has fewer than 2 children."
        )
    left = children[0].get_terminals()[0].name
    right = children[1].get_terminals()[0].name
    return (left, right)


def root_spanning_taxa(source: Path | str) -> tuple[str, str]:
    """Return two tips spanning the root of *source* (a path or Newick string).

    Used to anchor the root age when no species-pair calibration is supplied,
    which lets r8s produce a relative-time ultrametric tree.
    """
    tree = Phylo.read(StringIO(_read_newick_text(source)), "newick")
    return _root_spanning_taxa(tree)


# ===================================================================
#  4.  r8s control file + command
# ===================================================================


def _fmt_number(value: float) -> str:
    """Format a number without a spurious trailing ``.0`` (e.g. 94.0 -> ``94``)."""
    return f"{value:g}"


def _calibration_lines(calibrations: list[Calibration]) -> tuple[list[str], list[str]]:
    """Build the ``mrca`` lines and the ``fixage`` / ``constrain`` lines."""
    mrca_lines: list[str] = []
    age_lines: list[str] = []
    for cal in calibrations:
        if len(cal.taxa) != 2:
            raise ValueError(
                f"Calibration '{cal.name}' must name exactly two tip labels; "
                f"got {cal.taxa!r}"
            )
        if cal.age is None and cal.min_age is None and cal.max_age is None:
            raise ValueError(
                f"Calibration '{cal.name}' needs a fixed 'age' or a "
                "'min_age'/'max_age' window."
            )
        sp1, sp2 = cal.taxa
        mrca_lines.append(f"mrca {cal.name} {sp1} {sp2};")
        if cal.age is not None:
            age_lines.append(f"fixage taxon={cal.name} age={_fmt_number(cal.age)};")
        else:
            parts = [f"constrain taxon={cal.name}"]
            if cal.min_age is not None:
                parts.append(f"min_age={_fmt_number(cal.min_age)}")
            if cal.max_age is not None:
                parts.append(f"max_age={_fmt_number(cal.max_age)}")
            age_lines.append(" ".join(parts) + ";")
    return mrca_lines, age_lines


def _divtime_lines(
    method: str, algorithm: str, smoothing: float, cross_validate: bool
) -> list[str]:
    """Build the r8s smoothing / ``divtime`` command line(s).

    Default (``cross_validate=False``): a single penalized-likelihood fit at a
    fixed smoothing level. Cross-validation re-optimizes the whole tree once per
    terminal per smoothing value, so it scales as ~O(taxa) and is only practical
    for small trees — hence it is opt-in.
    """
    if cross_validate:
        return [
            f"divtime method={method} algorithm={algorithm} "
            "cvStart=0 cvInc=0.5 cvNum=8 crossv=yes;"
        ]
    lines: list[str] = []
    if method == "pl":
        lines.append(f"set smoothing={_fmt_number(smoothing)};")
    lines.append(f"divtime method={method} algorithm={algorithm};")
    return lines


def build_r8s_control(
    tree_newick: str,
    nsites: int,
    calibrations: Iterable[Calibration],
    *,
    method: str = "pl",
    algorithm: str = "tn",
    smoothing: float = 100.0,
    cross_validate: bool = False,
) -> str:
    """Build the NEXUS-style r8s control file text.

    Based on the CAFE-5 tutorial's ``prep_r8s.py``, generalised to support
    multiple calibrations, windowed constraints, and a configurable smoothing /
    cross-validation policy.

    Parameters
    ----------
    tree_newick : str
        A bare rooted Newick (typically from :func:`sanitize_tree_for_r8s`),
        ending in ``;``.
    nsites : int
        Number of alignment columns (see :func:`count_alignment_sites`).
    calibrations : iterable of Calibration
        At least one calibration.  At least one fixed age (or a fully bounded
        ``min_age``+``max_age`` window) is needed for r8s to set an absolute
        timescale; otherwise a warning is logged.
    method, algorithm : str
        r8s ``divtime`` options (defaults ``pl`` / ``tn``).
    smoothing : float
        Penalized-likelihood smoothing level used for a single fit (``set
        smoothing=…``). Ignored when ``cross_validate`` is True or ``method`` is
        not ``pl``.
    cross_validate : bool
        If True, cross-validate the smoothing parameter instead of a single fit.
        Accurate but scales as ~O(taxa) (the tutorial default) — only practical
        for small trees, so it is opt-in.

    Returns
    -------
    str
        The control-file text (newline-terminated).
    """
    if nsites <= 0:
        raise ValueError(f"nsites must be a positive integer, got {nsites}")

    calibrations = list(calibrations)
    if not calibrations:
        raise ValueError(
            "At least one calibration is required to make the tree ultrametric."
        )

    tree_line = tree_newick.strip()
    if not tree_line.endswith(";"):
        tree_line += ";"

    mrca_lines, age_lines = _calibration_lines(calibrations)

    pins_scale = any(
        cal.age is not None or (cal.min_age is not None and cal.max_age is not None)
        for cal in calibrations
    )
    if not pins_scale:
        logger.warning(
            "No calibration provides a fixed age or a bounded (min+max) window; "
            "r8s may be unable to set an absolute timescale."
        )

    divtime_lines = _divtime_lines(method, algorithm, smoothing, cross_validate)

    lines = [
        "#NEXUS",
        "begin trees;",
        f"tree {_TREE_NAME} = [&R] {tree_line}",
        "End;",
        "begin rates;",
        f"blformat nsites={nsites} lengths=persite ultrametric=no;",
        "collapse;",
        *mrca_lines,
        *age_lines,
        *divtime_lines,
        "describe plot=chronogram;",
        "describe plot=tree_description;",
        "end;",
        "",
    ]
    return "\n".join(lines)


def build_command(control_path: Path | str, tool_path: str = "r8s") -> list[str]:
    """Build the r8s command: ``r8s -b -f <control_path>``.

    ``-b`` runs in batch mode (suppresses the interactive banner).  r8s writes
    the dated tree to **stdout**, so the caller captures stdout and passes it to
    :func:`parse_ultrametric_tree`.
    """
    return [str(tool_path), "-b", "-f", str(control_path)]


# ===================================================================
#  5.  Parsing r8s output
# ===================================================================


def parse_ultrametric_tree(r8s_output: str | Path) -> str:
    """Extract the dated (ultrametric) Newick from r8s output.

    Scans from the end for the last ``tree <name> = <newick>;`` line — which the
    ``describe plot=tree_description`` command emits after the input-tree echo
    and the chronogram plot — and returns the Newick, stripping any leading
    ``[&R]`` / ``[&U]`` rooting comment.  This is robust to the exact column
    offset, unlike the tutorial's ``cut -c 16-``.

    Raises
    ------
    ValueError
        If no tree-description line is found (usually means r8s failed).
    """
    text = (
        r8s_output.read_text(encoding="utf-8")
        if isinstance(r8s_output, Path)
        else r8s_output
    )
    for line in reversed(text.splitlines()):
        match = _TREE_LINE_RE.search(line)
        if match:
            newick = match.group(1).strip()
            return _ROOTING_COMMENT_RE.sub("", newick)
    raise ValueError(
        "Could not find a tree-description line ('tree ... = ...;') in r8s output. "
        "r8s likely failed — check its stdout/stderr log."
    )


# ===================================================================
#  6.  Orchestration
# ===================================================================


def make_ultrametric(
    orthofinder_output_dir: Path | str | None,
    calibrations: Iterable[Calibration],
    *,
    out_tree: Path | str,
    work_dir: Path | str,
    input_tree: Path | str | None = None,
    nsites: int | None = None,
    r8s_path: str = "r8s",
    method: str = "pl",
    algorithm: str = "tn",
    smoothing: float = 100.0,
    cross_validate: bool = False,
    root_age: float = 1.0,
    dry_run: bool = False,
) -> dict:
    """Run the full OrthoFinder-tree → ultrametric-tree step via r8s.

    Steps: locate the species tree + alignment → derive ``nsites`` (unless
    overridden) → verify calibration taxa are tips → sanitize the tree → write
    the r8s control file → run r8s → parse and write the dated tree.

    All external-tool invocation goes through
    :func:`convgeno.utils.command_runner.run` (imported lazily, since that
    module is implemented in a later step of this plan).

    Parameters
    ----------
    orthofinder_output_dir : Path, str, or None
        OrthoFinder ``-o`` output dir, or a ``Results_*`` dir directly. May be
        ``None`` when ``input_tree`` is given together with an explicit
        ``nsites`` — then no OrthoFinder output is read at all.
    calibrations : iterable of Calibration
        Node calibrations. May be empty — then the root is anchored at
        ``root_age`` and r8s produces a *relative-time* ultrametric tree.
    out_tree : Path or str
        Where to write the ultrametric Newick.
    work_dir : Path or str
        Scratch dir for the control file and raw r8s output.
    input_tree : Path or str, optional
        Date THIS rooted Newick instead of discovering OrthoFinder's
        ``SpeciesTree_rooted.txt``. When given without ``nsites``, OrthoFinder's
        alignment is still read (from ``orthofinder_output_dir``) to derive it.
    nsites : int, optional
        Override the auto-derived alignment column count.
    r8s_path : str
        Path to (or name of) the r8s binary.
    method, algorithm, smoothing, cross_validate
        Passed to :func:`build_r8s_control`. Default is a single PL fit at
        ``smoothing``; set ``cross_validate=True`` only for small trees.
    root_age : float
        Root age used to anchor the tree when ``calibrations`` is empty
        (relative time; default 1.0).
    dry_run : bool
        Write the control file and log the command, but do not run r8s.

    Returns
    -------
    dict
        Stats: ``nsites``, ``n_calibrations``, ``tips``, ``r8s_ctl``,
        ``out_tree`` (and ``r8s_output`` on a real run, ``dry_run`` if dry).
    """
    # Imported here (not at module top) so convgeno.external.r8s stays importable
    # before convgeno.utils.command_runner is implemented.
    from convgeno.utils.command_runner import check_tool_available, run

    out_tree = Path(out_tree)
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    calibrations = list(calibrations)

    if input_tree is not None:
        species_tree_path = Path(input_tree)
        if not species_tree_path.is_file():
            raise FileNotFoundError(
                f"Input species tree not found: {species_tree_path}"
            )
        if nsites is None:
            # r8s scales per-site branch lengths by nsites; with a user tree but
            # no explicit nsites, fall back to OrthoFinder's alignment width.
            if orthofinder_output_dir is None:
                raise ValueError(
                    "nsites is required with input_tree when no OrthoFinder output "
                    "dir is given (pass nsites, or the OrthoFinder output dir)."
                )
            _, alignment_path = find_species_tree_files(orthofinder_output_dir)
            nsites = count_alignment_sites(alignment_path)
        else:
            logger.info(
                "Dating user-supplied tree %s with nsites=%d.",
                species_tree_path,
                nsites,
            )
    else:
        if orthofinder_output_dir is None:
            raise ValueError(
                "Either orthofinder_output_dir or input_tree must be given."
            )
        species_tree_path, alignment_path = find_species_tree_files(
            orthofinder_output_dir
        )
        if nsites is None:
            nsites = count_alignment_sites(alignment_path)
        else:
            logger.info(
                "Using caller-supplied nsites=%d (overriding alignment count).",
                nsites,
            )

    species_tree_obj = Phylo.read(str(species_tree_path), "newick")
    tip_names = {tip.name for tip in species_tree_obj.get_terminals()}

    relative_time = not calibrations
    if relative_time:
        anchor = _root_spanning_taxa(species_tree_obj)
        calibrations = [Calibration(name="root", taxa=anchor, age=root_age)]
        logger.warning(
            "No calibration supplied; producing a RELATIVE-time ultrametric tree "
            "anchored at root age=%s via taxa %s. Supply a calibration for "
            "absolute divergence times.",
            _fmt_number(root_age),
            anchor,
        )

    missing = {sp for cal in calibrations for sp in cal.taxa if sp not in tip_names}
    if missing:
        raise ValueError(
            f"Calibration taxa not found as tips in {species_tree_path.name}: "
            f"{sorted(missing)}. Tip labels are the OrthoFinder proteome basenames: "
            f"{sorted(tip_names)}"
        )

    sanitized = sanitize_tree_for_r8s(species_tree_path)
    control_text = build_r8s_control(
        sanitized,
        nsites,
        calibrations,
        method=method,
        algorithm=algorithm,
        smoothing=smoothing,
        cross_validate=cross_validate,
    )
    control_path = work_dir / "r8s_ctl_file.txt"
    control_path.write_text(control_text, encoding="utf-8")
    logger.info("Wrote r8s control file: %s", control_path)

    if not dry_run and not check_tool_available(r8s_path):
        raise FileNotFoundError(
            f"r8s not found or not executable: '{r8s_path}'. r8s is not on conda; "
            "build it with tools/r8s/install_r8s.sh (see tools/r8s/README.md), then "
            "put it on PATH or set its path in tool_paths.local.yaml."
        )

    cmd = build_command(control_path, tool_path=r8s_path)
    result = run(cmd, log_dir=work_dir, dry_run=dry_run)

    stats: dict = {
        "nsites": nsites,
        "n_calibrations": len(calibrations),
        "relative_time": relative_time,
        "tips": sorted(tip_names),
        "smoothing": smoothing,
        "cross_validate": cross_validate,
        "r8s_ctl": str(control_path),
        "out_tree": str(out_tree),
    }

    if dry_run:
        logger.info("Dry run: r8s was not executed.")
        stats["dry_run"] = True
        return stats

    r8s_stdout = result.stdout or ""
    raw_out = work_dir / "r8s_tmp.txt"
    raw_out.write_text(r8s_stdout, encoding="utf-8")

    newick = parse_ultrametric_tree(r8s_stdout)
    out_tree.parent.mkdir(parents=True, exist_ok=True)
    out_tree.write_text(newick + "\n", encoding="utf-8")
    logger.info("Wrote ultrametric species tree: %s", out_tree)

    stats["r8s_output"] = str(raw_out)
    return stats
