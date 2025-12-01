#!/usr/bin/env bash
set -euo pipefail

# Use vsearch to find the closest SILVA SSU match for each CoSMIC rRNA
# (rrna_uid from barrnap_rrna_sequences.fasta). CoSMIC then attaches these
# SILVA annotations to each ASV via the rrna_uid in the report.
#
# Usage:
#   ./rrna_vs_SILVA_vsearch.sh [RUN_OUTDIR] [SILVA_FASTA] [OUT_TSV]
# Defaults:
#   RUN_OUTDIR  = run_allMAGs-odSpoOffi2_id95-vsearch
#   SILVA_FASTA = silva/SILVA_138_SSURef_tax_silva.fasta
#   OUT_TSV     = SILVA_best_hits.tsv
#
# The OUT_TSV format is:
#   rrna_uid <TAB> silva_identity <TAB> silva_id <TAB> silva_taxonomy
# where silva_identity is a fraction (e.g. 0.987).

RUN_OUTDIR="${1:-run_allMAGs-odSpoOffi2_id95-vsearch}"
SILVA_FASTA="${2:-silva/SILVA_138_SSURef_tax_silva.fasta}"
OUT_TSV="${3:-SILVA_best_hits.tsv}"

RRNA_FASTA="${RUN_OUTDIR}/barrnap_rrna_sequences.fasta"

if [[ ! -f "$RRNA_FASTA" ]]; then
  echo "rRNA FASTA not found: $RRNA_FASTA" >&2
  echo "Make sure you ran the CoSMIC pipeline and that this OUTDIR is correct." >&2
  exit 1
fi

if [[ ! -f "$SILVA_FASTA" ]]; then
  echo "SILVA FASTA not found: $SILVA_FASTA" >&2
  exit 1
fi

if ! command -v vsearch >/dev/null 2>&1; then
  echo "vsearch not found on PATH. Ensure it is installed and available." >&2
  exit 1
fi

TMP_B6="$(mktemp)"
trap 'rm -f "$TMP_B6"' EXIT

MIN_ID="${VSEARCH_MIN_ID:-0.80}"
THREADS="${VSEARCH_THREADS:-8}"

echo "[SILVA] Running vsearch --usearch_global (min id=${MIN_ID}, threads=${THREADS})..." >&2
vsearch --usearch_global "$RRNA_FASTA" \
  --db "$SILVA_FASTA" \
  --id "$MIN_ID" \
  --strand both \
  --blast6out "$TMP_B6" \
  --threads "$THREADS"

python3 - "$TMP_B6" "$SILVA_FASTA" "$OUT_TSV" << 'PY'
import sys
from pathlib import Path

blast6_path = Path(sys.argv[1])
silva_fasta = Path(sys.argv[2])
out_path = Path(sys.argv[3])

if not blast6_path.exists():
    sys.exit(f"BLAST6 file not found: {blast6_path}")
if not silva_fasta.exists():
    sys.exit(f"SILVA FASTA not found: {silva_fasta}")

print(f"[SILVA] Parsing SILVA FASTA headers from {silva_fasta}...", file=sys.stderr)
silva_taxonomy = {}
with silva_fasta.open() as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        parts = line[1:].strip().split(" ", 1)
        accession = parts[0]
        taxonomy = parts[1] if len(parts) > 1 else ""
        silva_taxonomy[accession] = taxonomy

print(f"[SILVA] Parsing vsearch BLAST6 output from {blast6_path}...", file=sys.stderr)
best = {}
with blast6_path.open() as fh:
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip().split("\t")
        if len(parts) < 12:
            continue
        qid, sid, pid, aln_len, mism, gaps, qstart, qend, sstart, send, evalue, bitscore = parts[:12]
        try:
            pid_f = float(pid) / 100.0  # convert percent identity to fraction
            aln_i = int(aln_len)
            bits_f = float(bitscore)
        except ValueError:
            continue
        if aln_i <= 0:
            continue
        prev = best.get(qid)
        # Prefer higher identity; break ties with higher bitscore, then longer alignment
        if prev is None or pid_f > prev[0] or (
            pid_f == prev[0] and (bits_f > prev[2] or (bits_f == prev[2] and aln_i > prev[1]))
        ):
            best[qid] = (pid_f, aln_i, bits_f, sid)

print(f"[SILVA] Writing best hits to {out_path}...", file=sys.stderr)
with out_path.open("w") as out_fh:
    for qid, (pid_f, aln_i, bits_f, sid) in sorted(best.items(), key=lambda x: x[1][0], reverse=True):
        taxonomy = silva_taxonomy.get(sid, "")
        out_fh.write(f"{qid}\t{pid_f:.3f}\t{sid}\t{taxonomy}\n")

PY

echo "[SILVA] Done. Best-hit table written to: ${OUT_TSV}" >&2

