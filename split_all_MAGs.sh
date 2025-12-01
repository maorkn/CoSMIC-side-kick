#!/usr/bin/env bash
set -euo pipefail

# 1) split all_MAGs.fasta.gz into per-contig MAG_*.fasta
gzip -cd Data/all_MAGs.fasta.gz | python3 - << 'PY'
import sys
from pathlib import Path

out_dir = Path("Data")
out_dir.mkdir(exist_ok=True)

fh = None

def close():
    global fh
    if fh:
        fh.close()
        fh = None

for line in sys.stdin:
    if line.startswith(">"):
        close()
        header = line[1:].strip().split()[0]      # e.g. ENA|CAXXXX010000001|CAXXXX010000001.1
        mag_id = header.replace("|", "_")         # sanitize for filename
        out_path = out_dir / f"MAG_{mag_id}.fasta"
        fh = out_path.open("w")
    if fh:
        fh.write(line)

close()
PY

# 2) compress the new MAG FASTAs
for f in Data/MAG_*.fasta; do
  gzip "$f"
done
