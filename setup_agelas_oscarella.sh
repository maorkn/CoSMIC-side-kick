#!/bin/bash
# Setup script: extract zips and prepare data directories for Agelas and Oscarella runs
set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$BASEDIR"

echo "=== Extracting Agelas zip ==="
if [ ! -d "Agelas_data" ]; then
    mkdir -p Agelas_data
    unzip -o "Agelas and Oscarella/metagenoms Agelas.zip" -d Agelas_data/
    # Move files out of the subfolder
    if [ -d "Agelas_data/metagenoms Agelas" ]; then
        mv "Agelas_data/metagenoms Agelas"/* Agelas_data/ 2>/dev/null || true
        rmdir "Agelas_data/metagenoms Agelas" 2>/dev/null || true
    fi
    echo "Agelas extracted to Agelas_data/"
else
    echo "Agelas_data/ already exists, skipping extraction"
fi

echo ""
echo "=== Extracting Oscarella zip ==="
if [ ! -d "Oscarella_data" ]; then
    mkdir -p Oscarella_data
    unzip -o "Agelas and Oscarella/oscarella metagenome.zip" -d Oscarella_data/
    # Move files out of the subfolder
    if [ -d "Oscarella_data/oscarella metagenome" ]; then
        mv "Oscarella_data/oscarella metagenome"/* Oscarella_data/ 2>/dev/null || true
        rmdir "Oscarella_data/oscarella metagenome" 2>/dev/null || true
    fi
    echo "Oscarella extracted to Oscarella_data/"
else
    echo "Oscarella_data/ already exists, skipping extraction"
fi

# Create MAG-only subdirectories (only .fasta.gz and .fa.gz files, excluding metabarcoding fastas)
echo ""
echo "=== Creating MAG directories ==="

mkdir -p Agelas_data/MAGs
for f in Agelas_data/*.fasta.gz Agelas_data/*.fa.gz; do
    [ -f "$f" ] || continue
    bn=$(basename "$f")
    # Skip the metabarcoding / filtered sequences
    case "$bn" in
        *filtered*|*metagenome*) continue ;;
    esac
    ln -sf "$(realpath "$f")" "Agelas_data/MAGs/$bn" 2>/dev/null || cp "$f" "Agelas_data/MAGs/$bn"
done
echo "Agelas MAGs linked: $(ls Agelas_data/MAGs/ | wc -l) files"

mkdir -p Oscarella_data/MAGs
for f in Oscarella_data/*.fasta.gz Oscarella_data/*.fa.gz Oscarella_data/*.fasta; do
    [ -f "$f" ] || continue
    bn=$(basename "$f")
    # Skip the metabarcoding / filtered sequences and plain text files
    case "$bn" in
        *filtered*|*metagenome*|*.txt|*.json|*.tsv) continue ;;
    esac
    ln -sf "$(realpath "$f")" "Oscarella_data/MAGs/$bn" 2>/dev/null || cp "$f" "Oscarella_data/MAGs/$bn"
done
echo "Oscarella MAGs linked: $(ls Oscarella_data/MAGs/ | wc -l) files"

echo ""
echo "=== Listing metabarcoding FASTA files ==="
echo "Agelas full: $(ls -la Agelas_data/Agelas_filtered_sequences.fasta 2>/dev/null || echo 'NOT FOUND')"
echo "Agelas V3V4: $(ls -la Agelas_data/Agelas_filtered_V3V4_sequences.fasta 2>/dev/null || echo 'NOT FOUND')"
echo "Oscarella full: $(ls -la Oscarella_data/Oscarella_filtered_sequences.fasta 2>/dev/null || echo 'NOT FOUND')"
echo "Oscarella V3V4: $(ls -la Oscarella_data/Oscarella_filtered_V3V4_sequences.fasta 2>/dev/null || echo 'NOT FOUND')"

echo ""
echo "=== Setup complete ==="
