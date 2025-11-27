#!/usr/bin/env bash
set -euo pipefail

# Usage: ./summarize_GCA_vs_SILVA.sh [PAF_FILE] [SILVA_FASTA]
# Default PAF file is GCA_vs_SILVA.paf
# Default SILVA FASTA is silva/SILVA_138_SSURef_tax_silva.fasta
PAF_FILE="${1:-GCA_vs_SILVA.paf}"
SILVA_FASTA="${2:-silva/SILVA_138_SSURef_tax_silva.fasta}"

if [[ ! -f "$PAF_FILE" ]]; then
  echo "PAF file not found: $PAF_FILE" >&2
  exit 1
fi

if [[ ! -f "$SILVA_FASTA" ]]; then
  echo "SILVA FASTA file not found: $SILVA_FASTA" >&2
  exit 1
fi

python3 - "$PAF_FILE" "$SILVA_FASTA" << 'PY'
import sys
from pathlib import Path

paf_file = Path(sys.argv[1])
silva_fasta = Path(sys.argv[2])

# Parse SILVA FASTA headers to extract taxonomy
silva_taxonomy = {}
with silva_fasta.open() as fh:
    for line in fh:
        if line.startswith('>'):
            # Format: >ACCESSION TAXONOMY
            parts = line[1:].strip().split(' ', 1)
            accession = parts[0]
            taxonomy = parts[1] if len(parts) > 1 else ''
            silva_taxonomy[accession] = taxonomy

# Process PAF file
best = {}
with paf_file.open() as fh:
    for line in fh:
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.rstrip().split('\t')
        if len(parts) < 12:
            continue
        q, qlen, qs, qe, strand, t, tlen, ts, te, nm, aln, mapq = parts[:12]
        try:
            nm_i = int(nm)
            aln_i = int(aln)
        except ValueError:
            continue
        if aln_i <= 0:
            continue
        ident = nm_i / aln_i
        prev = best.get(q)
        if (not prev) or ident > prev[0]:
            best[q] = (ident, t)

# Output with taxonomy
for q, (ident, t) in sorted(best.items(), key=lambda x: x[1][0], reverse=True):
    taxonomy = silva_taxonomy.get(t, '')
    print(f"{q}\t{ident:.3f}\t{t}\t{taxonomy}")
PY
