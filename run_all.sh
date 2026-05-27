#!/bin/bash
# Master script: Run the CoSMIC pipeline for all 4 configurations
# Agelas V3V4, Agelas full, Oscarella V3V4, Oscarella full
set -euo pipefail

source ~/miniforge3/bin/activate cosmic_sidekick_env

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$BASEDIR"

echo "============================================="
echo "  CoSMIC Pipeline: Agelas & Oscarella Runs"
echo "============================================="
echo ""

# --- Run 1: Agelas V3V4 ---
echo "===== [1/4] Agelas V3V4 ====="
python3 cosmic_sidekick.py run \
  --config config_agelas_V3V4.yaml \
  --output-dir run_Agelas_V3V4
echo ""
echo "Generating report for Agelas V3V4..."
python3 cosmic_sidekick.py report \
  --config config_agelas_V3V4.yaml \
  --output-dir run_Agelas_V3V4
echo "===== Agelas V3V4 DONE ====="
echo ""

# --- Run 2: Agelas full ---
echo "===== [2/4] Agelas full-length 16S ====="
python3 cosmic_sidekick.py run \
  --config config_agelas_full.yaml \
  --output-dir run_Agelas_full
echo ""
echo "Generating report for Agelas full..."
python3 cosmic_sidekick.py report \
  --config config_agelas_full.yaml \
  --output-dir run_Agelas_full
echo "===== Agelas full DONE ====="
echo ""

# --- Run 3: Oscarella V3V4 ---
echo "===== [3/4] Oscarella V3V4 ====="
python3 cosmic_sidekick.py run \
  --config config_oscarella_V3V4.yaml \
  --output-dir run_Oscarella_V3V4
echo ""
echo "Generating report for Oscarella V3V4..."
python3 cosmic_sidekick.py report \
  --config config_oscarella_V3V4.yaml \
  --output-dir run_Oscarella_V3V4
echo "===== Oscarella V3V4 DONE ====="
echo ""

# --- Run 4: Oscarella full ---
echo "===== [4/4] Oscarella full-length 16S ====="
python3 cosmic_sidekick.py run \
  --config config_oscarella_full.yaml \
  --output-dir run_Oscarella_full
echo ""
echo "Generating report for Oscarella full..."
python3 cosmic_sidekick.py report \
  --config config_oscarella_full.yaml \
  --output-dir run_Oscarella_full
echo "===== Oscarella full DONE ====="
echo ""

echo "============================================="
echo "  All 4 runs complete!"
echo "============================================="
echo "Outputs:"
echo "  run_Agelas_V3V4/"
echo "  run_Agelas_full/"
echo "  run_Oscarella_V3V4/"
echo "  run_Oscarella_full/"
