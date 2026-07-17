#!/usr/bin/env bash
# Measure DIAMOND search throughput in bytes^2 per core-second.
#
# This is the per-core rate `throughput_const` used by the search-array
# walltime model in convgeno.slurm.multinode_generator.derive_search_walltime:
#
#     seconds = max_bucket_cost / (cores * throughput_const) * margin
#
# Cost model (identical to the generator):
#     cost(i, j)       = bytes(Species_i.fa) * bytes(Species_j.fa)
#     throughput_const = sum(cost) / sum(elapsed_seconds)
#
# Each sampled emitted DIAMOND command is run sequentially with ONE thread
# (-p 1); its output is redirected to a throwaway temp path so the real
# OrthoFinder WorkingDirectory search results are never modified.
#
# Usage:
#   bash throughput_calibration.sh SEARCH_COMMANDS WORK_DIR [N]
#
#   SEARCH_COMMANDS  the prepare job's <RUN>_search_commands.txt
#   WORK_DIR         the OrthoFinder WorkingDirectory (has Species*.fa)
#   N                number of commands to sample (default 100)
#
# Prereq: `diamond` on PATH (activate the convgeno env first).

set -uo pipefail

SEARCH_COMMANDS="${1:?Usage: throughput_calibration.sh SEARCH_COMMANDS WORK_DIR [N]}"
WORK_DIR="${2:?Usage: throughput_calibration.sh SEARCH_COMMANDS WORK_DIR [N]}"
N="${3:-100}"

command -v diamond >/dev/null 2>&1 || {
    echo "ERROR: diamond is not available on PATH (activate the convgeno env)." >&2
    exit 1
}
[ -f "$SEARCH_COMMANDS" ] || {
    echo "ERROR: search-command file not found: $SEARCH_COMMANDS" >&2
    exit 1
}
[ -d "$WORK_DIR" ] || {
    echo "ERROR: WorkingDirectory not found: $WORK_DIR" >&2
    exit 1
}

TOTAL=$(wc -l < "$SEARCH_COMMANDS")
if [ "$N" -gt "$TOTAL" ]; then
    N="$TOTAL"
fi

TEMP_DIR=$(mktemp -d)
RESULTS="$TEMP_DIR/results.tsv"
trap 'rm -rf "$TEMP_DIR"' EXIT

echo "Search commands:    $SEARCH_COMMANDS"
echo "WorkingDirectory:   $WORK_DIR"
echo "Available commands: $TOTAL"
echo "Sample size:        $N"
echo "DIAMOND:            $(command -v diamond)"
echo

shuf -n "$N" "$SEARCH_COMMANDS" > "$TEMP_DIR/sample.txt"

INDEX=0
SUCCESSFUL=0

while IFS= read -r COMMAND; do
    [ -n "$COMMAND" ] || continue
    INDEX=$((INDEX + 1))

    # Species pair from the version-stable Blast{i}_{j} output token.
    if [[ "$COMMAND" =~ Blast([0-9]+)_([0-9]+) ]]; then
        QUERY_ID="${BASH_REMATCH[1]}"
        DATABASE_ID="${BASH_REMATCH[2]}"
    else
        echo "[skip $INDEX] could not identify Blast{i}_{j}" >&2
        continue
    fi

    QUERY_FASTA="$WORK_DIR/Species${QUERY_ID}.fa"
    DATABASE_FASTA="$WORK_DIR/Species${DATABASE_ID}.fa"
    if [ ! -f "$QUERY_FASTA" ] || [ ! -f "$DATABASE_FASTA" ]; then
        echo "[skip $INDEX] missing Species FASTA file" >&2
        continue
    fi

    QUERY_BYTES=$(stat -c '%s' "$QUERY_FASTA")
    DATABASE_BYTES=$(stat -c '%s' "$DATABASE_FASTA")

    # Force single-thread and redirect output to scratch: strip any existing
    # -p/--threads and -o/--out, then append our own.
    RUN_COMMAND=$(printf '%s' "$COMMAND" |
        sed -E 's#[[:space:]](-p|--threads)([[:space:]]+|=)[0-9]+# #g' |
        sed -E 's#[[:space:]](-o|--out)([[:space:]]+|=)[^[:space:]]+# #g')
    OUTPUT_FILE="$TEMP_DIR/output_${INDEX}.txt.gz"
    RUN_COMMAND="$RUN_COMMAND -p 1 -o $(printf '%q' "$OUTPUT_FILE")"

    START_TIME=$(date +%s.%N)
    if ! ( cd "$WORK_DIR" && eval "$RUN_COMMAND" ) >/dev/null 2>&1; then
        echo "[skip $INDEX] DIAMOND command failed" >&2
        continue
    fi
    END_TIME=$(date +%s.%N)

    ELAPSED=$(awk -v s="$START_TIME" -v e="$END_TIME" 'BEGIN { printf "%.6f", e - s }')
    COST=$(python -c "print($QUERY_BYTES * $DATABASE_BYTES)")

    printf '%s\t%s\n' "$COST" "$ELAPSED" >> "$RESULTS"
    SUCCESSFUL=$((SUCCESSFUL + 1))
    printf '[%d/%d] Blast%s_%s cost=%s time=%ss\n' \
        "$INDEX" "$N" "$QUERY_ID" "$DATABASE_ID" "$COST" "$ELAPSED"
done < "$TEMP_DIR/sample.txt"

echo
echo "============================================================"
echo "Successful commands: $SUCCESSFUL"
echo "============================================================"

python - "$RESULTS" <<'PY'
from pathlib import Path
import statistics
import sys

results_path = Path(sys.argv[1])
if not results_path.exists() or results_path.stat().st_size == 0:
    raise SystemExit("ERROR: no successful calibration results.")

costs, times = [], []
for line in results_path.read_text().splitlines():
    cost_text, time_text = line.split("\t")
    costs.append(int(cost_text))
    times.append(float(time_text))

total_cost, total_time = sum(costs), sum(times)
aggregate = total_cost / total_time
per_command = sorted(c / t for c, t in zip(costs, times) if t > 0)

print(f"commands                    : {len(costs)}")
print(f"sum(cost) [bytes^2]         : {total_cost:.6e}")
print(f"sum(time) [seconds]         : {total_time:.6f}")
print()
print(f">>> throughput_const = {aggregate:.6e} bytes^2 per core-second")
print()
print(
    "per-command rate min/median/max: "
    f"{per_command[0]:.6e} / "
    f"{statistics.median(per_command):.6e} / "
    f"{per_command[-1]:.6e}"
)
PY
