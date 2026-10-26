#!/bin/bash
# EXPERIMENTAL. Resumable Sage re-search of the PXD052932 mzML runs listed in RUNS_TSV; see
# README.md for every setting and output.
#
# usage: driver.sh [--dry-run] WORK_DIR
#
# A background fetch loop streams each run's mzML (gzipped) into WORK_DIR/work, staying at most
# MAX_RUNS_ON_DISK runs ahead of the search loop and waiting while free disk is below
# MIN_FREE_GB. The search loop runs Sage on each run, keeps its rank-1 PSMs and appends a
# receipt to WORK_DIR/results/run_manifest.tsv, then deletes the mzML. A run with a receipt is
# never searched again, so the driver can be stopped and restarted at any time.
set -u

HERE=$(cd "$(dirname "$0")" && pwd)
DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then DRY_RUN=1; shift; fi
WORK_DIR=${1:-${WORK_DIR:-}}
if [ -z "$WORK_DIR" ]; then
  echo "usage: $0 [--dry-run] WORK_DIR   (or set WORK_DIR)" >&2
  exit 2
fi
WORK_DIR=$(cd "$WORK_DIR" 2>/dev/null && pwd) || { echo "WORK_DIR doesn't exist: ${1:-$WORK_DIR}" >&2; exit 2; }

RUNS_TSV=${RUNS_TSV:-$WORK_DIR/runs.tsv}
DB_FASTA=${DB_FASTA:-$WORK_DIR/db/dog_db.fasta}
SAGE=${SAGE:-sage}
PYTHON=${PYTHON:-python3}
SAGE_PARAMS=${SAGE_PARAMS:-$HERE/sage_params.json}
NCONN=${NCONN:-16}
SAGE_THREADS=${SAGE_THREADS:-6}
MIN_FREE_GB=${MIN_FREE_GB:-5}
MAX_RUNS_ON_DISK=${MAX_RUNS_ON_DISK:-2}
NICE=${NICE:-19}
RETRY_FAILED=${RETRY_FAILED:-0}

W=$WORK_DIR/work
R=$WORK_DIR/results
L=$WORK_DIR/logs
MAN=$R/run_manifest.tsv
LOG=$L/driver.log
PARAMS=$WORK_DIR/sage_params.resolved.json
LOCK=$WORK_DIR/.driver.lock
MIN_FREE_KB=$((MIN_FREE_GB * 1024 * 1024))

log() { local m; m="$(date '+%F %T') $*"; echo "$m"; [ "$DRY_RUN" = 1 ] || echo "$m" >> "$LOG"; }
die() { echo "error: $*" >&2; exit 1; }
is_done() { [ -f "$MAN" ] && cut -f1 "$MAN" | grep -qxF "$1"; }
has_marker() { [ -f "$W/$1.fetch_failed" ] || [ -f "$W/$1.sage_failed" ] || [ -f "$W/$1.post_failed" ]; }
free_kb() { df -Pk "$WORK_DIR" | awk 'NR==2 {print $4}'; }
n_mzml() { find "$W" -maxdepth 1 \( -name '*.mzML.gz' -o -name '*.mzML.gz.part' \) | wc -l | tr -d ' '; }

# ---------------------------------------------------------------- checks
[ -f "$RUNS_TSV" ] || die "RUNS_TSV not found: $RUNS_TSV (make it with list_runs.py)"
[ -f "$DB_FASTA" ] || die "DB_FASTA not found: $DB_FASTA (make it with build_db.py)"
[ -f "$SAGE_PARAMS" ] || die "SAGE_PARAMS not found: $SAGE_PARAMS"
command -v "$SAGE" > /dev/null || die "Sage not found: $SAGE (set SAGE to the sage binary)"
"$PYTHON" -I -c 'import pandas, pyarrow, inflate64' 2> /dev/null \
  || die "$PYTHON needs pandas, pyarrow and inflate64"
for v in NCONN SAGE_THREADS MIN_FREE_GB MAX_RUNS_ON_DISK NICE; do
  case ${!v} in ''|*[!0-9]*) die "$v must be a non-negative integer, not '${!v}'";; esac
done
[ "$MAX_RUNS_ON_DISK" -ge 1 ] || die "MAX_RUNS_ON_DISK must be at least 1"

# Sage parameters with the database path filled in, plus the database's SHA-256. A WORK_DIR
# holds one search: restarting with other parameters or another database (by content, not
# path) is refused once a run has finished.
resolved=$("$PYTHON" -I - "$SAGE_PARAMS" "$DB_FASTA" <<'PY'
import hashlib, json, os, sys
params = json.load(open(sys.argv[1]))
params["database"]["fasta"] = os.path.abspath(sys.argv[2])
h = hashlib.sha256()
with open(sys.argv[2], "rb") as f:
    for block in iter(lambda: f.read(1 << 20), b""):
        h.update(block)
print(json.dumps({"sage_params": params, "db_fasta_sha256": h.hexdigest()}, indent=1, sort_keys=True))
PY
) || die "couldn't resolve Sage parameters"
if [ -f "$WORK_DIR/search_config.json" ] && [ -s "$MAN" ]; then
  "$PYTHON" -I -c '
import json, sys
old, new = json.load(open(sys.argv[1])), json.loads(sys.argv[2])
for config in (old, new):
    config["sage_params"]["database"].pop("fasta")
sys.exit(old != new)' "$WORK_DIR/search_config.json" "$resolved" \
    || die "WORK_DIR was searched with other Sage parameters or another database (see $WORK_DIR/search_config.json); use a new WORK_DIR"
fi

n_listed=0; n_done=0; n_failed=0; n_pending=0; pending_bytes=0; first=
while IFS=$'\t' read -r run grp url member usize csize; do
  [ -n "$run" ] || continue
  n_listed=$((n_listed + 1))
  if is_done "$run"; then n_done=$((n_done + 1))
  elif has_marker "$run" && [ "$RETRY_FAILED" != 1 ]; then n_failed=$((n_failed + 1))
  else
    n_pending=$((n_pending + 1)); pending_bytes=$((pending_bytes + csize))
    [ -n "$first" ] || first="$run $url $member"
  fi
done < "$RUNS_TSV"

if [ "$DRY_RUN" = 1 ]; then
  cat <<EOF
dry run: nothing is downloaded, searched, written or deleted.
WORK_DIR          $WORK_DIR
RUNS_TSV          $RUNS_TSV ($n_listed runs: $n_done done, $n_failed with failure markers, $n_pending to search)
DB_FASTA          $DB_FASTA
SAGE              $(command -v "$SAGE") ($("$SAGE" --version 2>/dev/null))
SAGE_PARAMS       $SAGE_PARAMS
PYTHON            $(command -v "$PYTHON")
NCONN=$NCONN SAGE_THREADS=$SAGE_THREADS MIN_FREE_GB=$MIN_FREE_GB MAX_RUNS_ON_DISK=$MAX_RUNS_ON_DISK NICE=$NICE RETRY_FAILED=$RETRY_FAILED
to download       $(awk -v b="$pending_bytes" 'BEGIN {printf "%.1f GB compressed, about %.1f h at 5 MB/s", b / 1e9, b / 5e6 / 3600}')
free disk         $(($(free_kb) / 1024 / 1024)) GB (fetching pauses below $MIN_FREE_GB GB)
EOF
  if [ -n "$first" ]; then
    read -r run url member <<< "$first"
    echo "first run         $run"
    echo "  fetch: nice -n $NICE $PYTHON -I $HERE/fetch_member.py $url $member $W/$run.mzML.gz $NCONN"
    echo "  search: RAYON_NUM_THREADS=$SAGE_THREADS nice -n $NICE $SAGE --disable-telemetry-i-dont-want-to-improve-sage -o $W/${run}_sage $PARAMS $W/$run.mzML.gz"
    echo "  keep: $PYTHON -I $HERE/post_run.py $run ... $R/psms/$run.parquet $MAN"
  fi
  exit 0
fi

mkdir -p "$W" "$R/psms" "$R/runs" "$L"
mkdir "$LOCK" 2> /dev/null || die "another driver is using $WORK_DIR (remove $LOCK if none is)"
FPID=
cleanup() {
  if [ -n "$FPID" ]; then pkill -P "$FPID" 2> /dev/null; kill "$FPID" 2> /dev/null; fi
  rmdir "$LOCK" 2> /dev/null
}
trap cleanup EXIT
trap 'log "INTERRUPTED"; exit 130' INT TERM

echo "$resolved" > "$WORK_DIR/search_config.json"
"$PYTHON" -I -c 'import json, sys; json.dump(json.load(open(sys.argv[1]))["sage_params"], open(sys.argv[2], "w"), indent=1)' \
  "$WORK_DIR/search_config.json" "$PARAMS"

# Leftovers of an interrupted session: partial downloads, and receipts without their mzML.
rm -f "$W"/*.part "$W"/*.part.bad "$W"/*.receipt.tmp
for receipt in "$W"/*.receipt.json; do
  [ -e "$receipt" ] || continue
  run=$(basename "$receipt" .receipt.json)
  [ -f "$W/$run.mzML.gz" ] || rm -f "$receipt"
done
if [ "$RETRY_FAILED" = 1 ]; then
  for marker in "$W"/*.fetch_failed "$W"/*.sage_failed "$W"/*.post_failed; do
    [ -e "$marker" ] || continue
    run=${marker##*/}; run=${run%.*}
    if [ -d "$W/${run}_sage" ]; then rm -f "$W/${run}_sage"/*; rmdir "$W/${run}_sage"; fi
    rm -f "$marker"
  done
fi

fetch_loop() {
  while IFS=$'\t' read -r run grp url member usize csize; do
    [ -n "$run" ] || continue
    is_done "$run" && continue
    has_marker "$run" && continue
    [ -f "$W/$run.receipt.json" ] && [ -f "$W/$run.mzML.gz" ] && continue
    while [ "$(n_mzml)" -ge "$MAX_RUNS_ON_DISK" ]; do sleep 10; done
    waited=0
    while [ "$(free_kb)" -lt "$MIN_FREE_KB" ]; do
      if [ $((waited % 900)) -eq 0 ]; then
        log "LOW_DISK free=$(($(free_kb) / 1024))MB, waiting before $run (waited ${waited}s)"
      fi
      sleep 60; waited=$((waited + 60))
    done
    ok=0
    for attempt in 1 2 3; do
      if nice -n "$NICE" "$PYTHON" -I "$HERE/fetch_member.py" "$url" "$member" "$W/$run.mzML.gz" "$NCONN" \
          > "$W/$run.receipt.tmp" 2>> "$L/fetch_$run.log"; then
        mv "$W/$run.receipt.tmp" "$W/$run.receipt.json"; ok=1; break
      fi
      log "FETCH_RETRY $run attempt $attempt"
      rm -f "$W/$run.mzML.gz.part" "$W/$run.mzML.gz.part.bad" "$W/$run.receipt.tmp"
      sleep 60
    done
    if [ $ok = 1 ]; then
      log "FETCHED $run $(tr -d '\n' < "$W/$run.receipt.json" | cut -c1-400)"
    else
      log "ERROR fetch failed $run (see logs/fetch_$run.log)"; touch "$W/$run.fetch_failed"
    fi
  done < "$RUNS_TSV"
  log "FETCH_LOOP_END"
}

search_loop() {
  i=0
  while IFS=$'\t' read -r run grp url member usize csize; do
    [ -n "$run" ] || continue
    i=$((i + 1))
    is_done "$run" && continue
    has_marker "$run" && continue
    until [ -f "$W/$run.receipt.json" ] || [ -f "$W/$run.fetch_failed" ]; do
      if ! kill -0 "$FPID" 2> /dev/null && [ ! -f "$W/$run.receipt.json" ]; then
        log "ERROR fetch loop ended before $run was fetched"; return 1
      fi
      sleep 10
    done
    if [ -f "$W/$run.fetch_failed" ]; then log "SKIP $run (fetch failed)"; continue; fi
    out=$W/${run}_sage
    t0=$(date +%s)
    if RAYON_NUM_THREADS=$SAGE_THREADS nice -n "$NICE" "$SAGE" --disable-telemetry-i-dont-want-to-improve-sage \
        -o "$out" "$PARAMS" "$W/$run.mzML.gz" > "$L/sage_$run.log" 2>&1; then
      t1=$(date +%s)
      if summary=$(nice -n "$NICE" "$PYTHON" -I "$HERE/post_run.py" "$run" "$out/results.sage.tsv" \
          "$W/$run.receipt.json" $((t1 - t0)) "$R/psms/$run.parquet" "$MAN" 2>> "$L/post_$run.log"); then
        cp "$out/results.json" "$R/runs/$run.sage_results.json"
        cp "$W/$run.receipt.json" "$R/runs/$run.receipt.json"
        rm -f "$W/$run.mzML.gz" "$W/$run.receipt.json" "$out"/*
        rmdir "$out"
        log "DONE $i $run sage_s=$((t1 - t0)) free=$(($(free_kb) / 1024))MB $(echo "$summary" | cut -c1-600)"
      else
        log "ERROR post_run $run (see logs/post_$run.log)"
        rm -f "$W/$run.mzML.gz" "$W/$run.receipt.json"; touch "$W/$run.post_failed"
      fi
    else
      log "ERROR sage $run (see logs/sage_$run.log)"
      rm -f "$W/$run.mzML.gz" "$W/$run.receipt.json"; touch "$W/$run.sage_failed"
    fi
  done < "$RUNS_TSV"
  log "SEARCH_LOOP_END"
}

log "START pid=$$ runs=$n_listed done=$n_done failed=$n_failed to_search=$n_pending nconn=$NCONN threads=$SAGE_THREADS"
fetch_loop &
FPID=$!
search_status=0
fetch_status=0
search_loop || search_status=$?
wait "$FPID" || fetch_status=$?
FPID=
n_marked=$(find "$W" -maxdepth 1 \( -name '*.fetch_failed' -o -name '*.sage_failed' -o -name '*.post_failed' \) | wc -l | tr -d ' ')
n_rows=0
[ -f "$MAN" ] && n_rows=$(($(wc -l < "$MAN") - 1))
if [ "$search_status" -ne 0 ] || [ "$fetch_status" -ne 0 ] || [ "$n_marked" -ne 0 ]; then
  log "INCOMPLETE searched=$n_rows failed=$n_marked search_status=$search_status fetch_status=$fetch_status"
  exit 1
fi
log "ALL_DONE searched=$n_rows failed=0"
