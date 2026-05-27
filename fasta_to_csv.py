#!/usr/bin/env python3
"""Convert a metabarcoding FASTA file to CSV format for cosmic_sidekick.

Each sequence gets ID from the FASTA header and a dummy uniform abundance.
"""
import sys
from pathlib import Path

def fasta_to_csv(fasta_path: str, csv_path: str):
    records = []
    header = None
    seq_parts = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                # Use entire header (without '>') as ID
                header = line[1:].split()[0]  # first word as ID
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts)))

    n = len(records)
    print(f"  {fasta_path}: {n} sequences")

    # Write CSV with dummy abundance column
    with open(csv_path, "w") as out:
        out.write("id,sequence,dummy_abundance\n")
        for seq_id, seq in records:
            # Give each sequence equal abundance (1/n)
            out.write(f"{seq_id},{seq},{1.0/n:.6f}\n")

    print(f"  -> {csv_path}")

if __name__ == "__main__":
    conversions = [
        ("Agelas_data/Agelas_filtered_sequences.fasta",
         "Agelas_data/Agelas_full_metabarcoding.csv"),
        ("Agelas_data/Agelas_filtered_V3V4_sequences.fasta",
         "Agelas_data/Agelas_V3V4_metabarcoding.csv"),
        ("Oscarella_data/Oscarella_filtered_sequences.fasta",
         "Oscarella_data/Oscarella_full_metabarcoding.csv"),
        ("Oscarella_data/Oscarella_filtered_V3V4_sequences.fasta",
         "Oscarella_data/Oscarella_V3V4_metabarcoding.csv"),
    ]

    for fasta, csv in conversions:
        if Path(fasta).exists():
            fasta_to_csv(fasta, csv)
        else:
            print(f"  WARNING: {fasta} not found, skipping")
