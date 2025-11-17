# CoSMIC Sidekick: MAG–Metabarcoding Linking Pipeline

This repository contains a small Python pipeline that:

- scans a directory of MAG assemblies (PacBio long-read, one MAG per fasta(.gz) file),
- uses **Barrnap** to extract 16S and 18S rRNA sequences from each MAG,
- creates a CSV mapping rRNA sequences back to their MAG,
- takes a metabarcoding CSV (full-length 16S/18S with per-sample relative abundances),
- maps those metabarcoding sequences to MAG rRNA sequences using a sequence-identity threshold,
- and annotates only the MAGs that have at least one metabarcoding hit using **Prokka**.

The main entrypoint is `cosmic_sidekick.py`.

## 1. Dependencies

### Python packages

Install Python dependencies (tested with Python 3.12):

```bash
pip install -r requirements.txt
```

### External tools

The script expects the following tools to be available on your `PATH`:

- **Barrnap** – for rRNA prediction (16S/18S).
- **Prokka** – for MAG annotation (simpler alternative to Trinotate).

Example installation options (pick what matches your environment):

- Using conda (recommended):

  ```bash
  conda install -c bioconda barrnap prokka
  ```

- Or using your distro’s package manager if available (names may vary):

  ```bash
  sudo apt-get update
  sudo apt-get install barrnap prokka
  ```

If `barrnap` or `prokka` are not found, the script will exit (for Barrnap) or skip annotation (for Prokka) with a clear message.

## 2. Configuration

The pipeline is configured via `config.yaml`. Key fields:

- `mags_dir`: directory with MAG assemblies (default: `Data`).
- `metabarcoding_csv`: path to your metabarcoding CSV.
- `metabarcoding_id_column`: column with unique IDs per metabarcoding sequence (e.g. ASV/OTU).
- `metabarcoding_sequence_column`: column with nucleotide sequences.
- `metabarcoding_abundance_columns`:
  - if `null`, the script auto-detects all numeric columns (except ID and sequence) as abundance columns;
  - or provide an explicit list of column names.
- `identity_threshold`: sequence identity threshold for mapping (e.g. `0.97` for 97%).
- `barrnap_kingdoms`: kingdoms to scan with Barrnap (default: `bac`, `arc`, `euk`).
- `annotation_tool`: currently implemented: `prokka`.
- `annotation_output_dir`: where MAG annotations are written (default: `Annotation/`).

Adjust `config.yaml` to match your metabarcoding CSV (especially the ID, sequence, and abundance columns) before running.

## 3. Running the pipeline

Basic usage:

```bash
python cosmic_sidekick.py run --mags-dir Data --metabarcoding your_metabarcoding.csv
```

or, relying on `config.yaml`:

```bash
python cosmic_sidekick.py run --config config.yaml
```

Command-line arguments override values in `config.yaml`.

## 4. Outputs

Running the full pipeline produces:

- `barrnap_rrna_mapping.csv` – one row per 16S/18S rRNA feature, with:
  - MAG ID (from filename),
  - contig/sequence header,
  - rRNA ID (from Barrnap),
  - rRNA type (`16S` / `18S`),
  - coordinates, strand, length,
  - kingdom, product,
  - extracted nucleotide sequence.
- `barrnap_rrna_sequences.fasta` – FASTA of all extracted rRNA sequences with stable IDs.
- `metabarcoding_to_MAG_mapping.csv` – rows mapping metabarcoding sequences to MAG rRNA hits:
  - metabarcoding ID and sequence,
  - MAG ID and rRNA info,
  - percent identity,
  - all per-sample abundance columns carried through.
- `Annotation/<MAG_ID>/` – Prokka outputs for each MAG that has at least one metabarcoding hit.

## 5. Typical workflow

1. Place your MAG assemblies in `Data/` as `.fasta` or `.fasta.gz` (one MAG per file).
2. Prepare your metabarcoding CSV (full-length 16S/18S, relative abundances per sample).
3. Edit `config.yaml` so `metabarcoding_csv`, `metabarcoding_id_column`,
   `metabarcoding_sequence_column`, and (optionally) `metabarcoding_abundance_columns`
   match your file.
4. Ensure `barrnap` and `prokka` are installed and on your `PATH`.
5. Run:

   ```bash
   python cosmic_sidekick.py run --config config.yaml
   ```

6. Inspect:
   - `barrnap_rrna_mapping.csv` / `barrnap_rrna_sequences.fasta` for rRNA extraction,
   - `metabarcoding_to_MAG_mapping.csv` for mapping results,
   - `Annotation/` for MAG annotations.

## 6. Notes and extensions

- Mapping uses a simple global identity measure (including a fast path for equal-length sequences, and a more general alignment for length mismatches). If your dataset becomes large, we can switch to `vsearch` or `blastn` for faster similarity searches.
- Currently, only Prokka-based annotation is implemented, but the code is structured so another tool (including Trinotate) could be plugged in if needed.

