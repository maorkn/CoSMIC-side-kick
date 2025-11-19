"""
Optional script to enrich the CoSMIC LLM report with additional
functional context from external tools.

It can:
1. Run default functional annotation for MAGs using:
   - DRAM (genome & metabolism summaries)
   - METABOLIC (metabolic marker summary)
   - eggNOG-mapper (KOs / GO terms from protein FASTAs)
2. Integrate these tool outputs (plus optional HUMAnN) into an enriched
   markdown report suitable as LLM context.

Usage (basic enrichment only, assuming tools already ran elsewhere):

    source .venv/bin/activate
    python Richer_report.py \
        --base-report cosmic_llm_report.md \
        --output cosmic_llm_rich_report.md

Run tools and enrich in one go (paths are examples; adjust to your setup):

    python Richer_report.py \
        --base-report cosmic_llm_report.md \
        --mags-dir Data \
        --annotation-dir Annotation \
        --run-dram \
        --run-metabolic \
        --run-eggnog \
        --threads 8 \
        --output cosmic_llm_rich_report.md

Notes:
- This script assumes DRAM (DRAM.py), METABOLIC-G, and eggNOG-mapper
  (emapper.py) are already installed and configured (including their
  databases) and available on your PATH. It does not install databases.
- HUMAnN and other tools are still expected to be run externally.
"""

import argparse
import io
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd
from shutil import which


def read_text_lines(path: Path) -> List[str]:
    if not path.exists():
        raise SystemExit(f"Base report not found: {path}")
    with path.open() as fh:
        return fh.read().splitlines()


def find_executable(candidates: Sequence[str]) -> Optional[str]:
    for name in candidates:
        exe = which(name)
        if exe is not None:
            return exe
    return None


def run_command(cmd: Sequence[str], cwd: Optional[Path] = None) -> bool:
    print(f"[Richer_report] Running: {' '.join(cmd)}", file=sys.stderr)
    try:
        subprocess.run(
            cmd,
            check=True,
            cwd=str(cwd) if cwd is not None else None,
        )
        return True
    except subprocess.CalledProcessError as e:
        print(
            f"[Richer_report] Command failed with exit code {e.returncode}: "
            f"{' '.join(cmd)}",
            file=sys.stderr,
        )
        return False


def detect_id_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower_map:
            return lower_map[cand.lower()]
    return None


def summarize_dram_genome_summary(path: Path) -> List[str]:
    lines: List[str] = []
    if not path or not path.exists():
        return lines

    header: Optional[str] = None
    data_lines: List[str] = []
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").strip()
                continue
            data_lines.append(line)

    if header is None:
        return lines

    buffer = io.StringIO(header + "\n" + "".join(data_lines))
    df = pd.read_csv(buffer, sep="\t")
    id_col = detect_id_column(df, ["genome", "Genome", "bin", "Bin", "MAG", "mag_id"])
    if id_col is None:
        id_col = df.columns[0]

    lines.append("## Additional Functional Context: DRAM Genome Summary")

    for _, row in df.iterrows():
        genome_id = row[id_col]
        lines.append(f"### DRAM genome summary for {genome_id}")

        # Prefer some commonly useful columns if present
        preferred_cols = [
            "Completeness",
            "Contamination",
            "size_Mbp",
            "coding_density",
            "tRNA",
            "rRNA",
        ]
        printed_any = False
        for col in preferred_cols:
            if col in df.columns:
                val = row[col]
                lines.append(f"- {col}: {val}")
                printed_any = True

        # As a fallback, show a few other numeric columns
        if not printed_any:
            numeric_cols = [
                c
                for c in df.columns
                if pd.api.types.is_numeric_dtype(df[c]) and c != id_col
            ]
            for col in numeric_cols[:10]:
                val = row[col]
                lines.append(f"- {col}: {val}")

        lines.append("")

    return lines


def summarize_dram_metabolism(path: Path) -> List[str]:
    lines: List[str] = []
    if not path or not path.exists():
        return lines

    df = pd.read_csv(path, sep="\t")
    id_col = detect_id_column(df, ["genome", "Genome", "bin", "Bin", "MAG", "mag_id"])
    if id_col is None:
        id_col = df.columns[0]

    lines.append("## Additional Functional Context: DRAM Metabolism Summary")

    for _, row in df.iterrows():
        genome_id = row[id_col]
        lines.append(f"### DRAM metabolism summary for {genome_id}")

        # Treat non-zero / non-empty values in remaining columns as "present"
        functional_cols = [c for c in df.columns if c != id_col]
        present_features: List[str] = []
        for col in functional_cols:
            val = row[col]
            if pd.isna(val):
                continue
            if isinstance(val, (int, float)) and val != 0:
                present_features.append(col)
            elif isinstance(val, str) and val.strip() not in {"0", "-", "NA", "False", "false"}:
                present_features.append(col)

        if present_features:
            lines.append("- Notable metabolic features present (by DRAM):")
            for feat in present_features[:30]:
                lines.append(f"  - {feat}")
        else:
            lines.append("- No non-zero metabolic features detected in this summary.")

        lines.append("")

    return lines


def summarize_metabolic_summary(path: Path) -> List[str]:
    lines: List[str] = []
    if not path or not path.exists():
        return lines

    df = pd.read_csv(path, sep="\t")
    id_col = detect_id_column(df, ["Genome", "genome", "MAG", "mag_id"])
    if id_col is None:
        id_col = df.columns[0]

    lines.append("## Additional Functional Context: METABOLIC Summary")

    for _, row in df.iterrows():
        genome_id = row[id_col]
        lines.append(f"### METABOLIC summary for {genome_id}")

        functional_cols = [c for c in df.columns if c != id_col]
        present_features: List[str] = []
        for col in functional_cols:
            val = row[col]
            if pd.isna(val):
                continue
            if isinstance(val, (int, float)) and val != 0:
                present_features.append(f"{col} (score={val})")
            elif isinstance(val, str) and val.strip() not in {"0", "-", "NA", "False", "false"}:
                present_features.append(f"{col}: {val}")

        if present_features:
            lines.append("- Selected metabolic markers (METABOLIC):")
            for feat in present_features[:30]:
                lines.append(f"  - {feat}")
        else:
            lines.append("- No non-zero metabolic markers detected in this summary.")

        lines.append("")

    return lines


def summarize_humann_pathways(path: Path) -> List[str]:
    lines: List[str] = []
    if not path or not path.exists():
        return lines

    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return lines

    # HUMAnN pathway tables typically: first column is pathway id/name
    pathway_col = df.columns[0]
    sample_cols = [c for c in df.columns[1:] if pd.api.types.is_numeric_dtype(df[c])]

    if not sample_cols:
        return lines

    df["total_abundance"] = df[sample_cols].sum(axis=1)
    top = df.sort_values("total_abundance", ascending=False).head(20)

    lines.append("## Additional Functional Context: HUMAnN Pathway Abundance")
    lines.append(
        "Top pathways across all samples (by total HUMAnN abundance):"
    )
    for _, row in top.iterrows():
        pid = row[pathway_col]
        total = row["total_abundance"]
        lines.append(f"- {pid}: total_abundance={total:.4f}")

    lines.append("")
    return lines


def summarize_eggnog_annotations(path: Path) -> List[str]:
    """
    Summarize eggNOG-mapper annotations from a combined TSV with columns like:
    'query', 'KEGG_ko', 'GOs', ...
    """
    lines: List[str] = []
    if not path or not path.exists():
        return lines

    header: Optional[str] = None
    data_lines: List[str] = []
    with path.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").strip()
                continue
            data_lines.append(line)

    if not header:
        return lines

    buffer = io.StringIO(header + "\n" + "".join(data_lines))
    df = pd.read_csv(buffer, sep="\t")

    ko_col = detect_id_column(df, ["KEGG_ko", "kegg_ko", "KO"])
    go_col = detect_id_column(df, ["GOs", "gos", "GO_terms"])

    lines.append("## Additional Functional Context: eggNOG-mapper")

    if ko_col and ko_col in df.columns:
        kos: Dict[str, int] = {}
        for val in df[ko_col].dropna():
            if not isinstance(val, str):
                continue
            for ko in [x.strip() for x in val.replace(";", ",").split(",")]:
                if not ko or ko in {"-", "NA"}:
                    continue
                kos[ko] = kos.get(ko, 0) + 1

        if kos:
            lines.append("- Top KEGG KOs (by annotation count):")
            for ko, count in sorted(kos.items(), key=lambda kv: kv[1], reverse=True)[:30]:
                lines.append(f"  - {ko} (n={count})")
        else:
            lines.append("- No KEGG KOs detected in eggNOG annotations.")

    if go_col and go_col in df.columns:
        gos: Dict[str, int] = {}
        for val in df[go_col].dropna():
            if not isinstance(val, str):
                continue
            for go in [x.strip() for x in val.replace(";", ",").split(",")]:
                if not go or go in {"-", "NA"}:
                    continue
                gos[go] = gos.get(go, 0) + 1

        if gos:
            lines.append("- Top GO terms (by annotation count):")
            for go, count in sorted(gos.items(), key=lambda kv: kv[1], reverse=True)[:30]:
                lines.append(f"  - {go} (n={count})")
        else:
            lines.append("- No GO terms detected in eggNOG annotations.")

    lines.append("")
    return lines


def run_dram(
    mags_dir: Path,
    out_dir: Path,
    threads: int,
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Run DRAM annotate + distill on MAG FASTA files.

    Returns (genome_summary.tsv path, metabolism_summary.tsv path), or (None, None)
    on failure.
    """
    exe = find_executable(["DRAM.py", "dram.py", "DRAM"])
    if exe is None:
        print(
            "[Richer_report] DRAM not found on PATH (DRAM.py/dram.py). "
            "Skipping DRAM run.",
            file=sys.stderr,
        )
        return None, None

    mags_dir = mags_dir.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    annotate_cmd = [
        exe,
        "annotate",
        "-i",
        str(mags_dir),
        "-o",
        str(out_dir),
        "--threads",
        str(threads),
    ]

    if not run_command(annotate_cmd):
        return None, None

    annotations_tsv = out_dir / "annotations.tsv"
    if not annotations_tsv.exists():
        print(
            f"[Richer_report] DRAM annotations.tsv not found in {out_dir}; "
            "cannot run distill.",
            file=sys.stderr,
        )
        return None, None

    distill_dir = out_dir / "distill"
    distill_cmd = [
        exe,
        "distill",
        "-i",
        str(annotations_tsv),
        "-o",
        str(distill_dir),
    ]
    if not run_command(distill_cmd):
        return None, None

    genome_summary = distill_dir / "genome_summary.tsv"
    metabolism_summary = distill_dir / "metabolism_summary.tsv"

    if not genome_summary.exists():
        genome_summary = None
    if not metabolism_summary.exists():
        metabolism_summary = None

    return genome_summary, metabolism_summary


def run_metabolic(
    mags_dir: Path,
    out_dir: Path,
    threads: int,
) -> Optional[Path]:
    """
    Run METABOLIC-G on MAG FASTA files.

    Returns a summary TSV path if found, otherwise None.
    """
    exe = find_executable(["METABOLIC-G", "METABOLIC-G.pl"])
    if exe is None:
        print(
            "[Richer_report] METABOLIC-G not found on PATH. "
            "Skipping METABOLIC run.",
            file=sys.stderr,
        )
        return None

    mags_dir = mags_dir.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    genome_list = out_dir / "METABOLIC_genomes.txt"
    mag_files = sorted(
        p
        for p in mags_dir.iterdir()
        if p.suffix in {".fa", ".fna", ".fasta", ".gz"}
    )
    if not mag_files:
        print(
            f"[Richer_report] No MAG FASTA files found in {mags_dir}; "
            "skipping METABOLIC run.",
            file=sys.stderr,
        )
        return None

    with genome_list.open("w") as fh:
        for p in mag_files:
            fh.write(str(p.resolve()) + "\n")

    cmd = [
        exe,
        "-in-gn",
        str(genome_list),
        "-o",
        str(out_dir),
        "-t",
        str(threads),
    ]

    if not run_command(cmd):
        return None

    # Try to find a summary TSV; METABOLIC naming conventions may vary slightly.
    candidate = out_dir / "METABOLIC_result_summary.tsv"
    if candidate.exists():
        return candidate

    # Fallback: first TSV file with 'summary' in the name.
    for p in out_dir.glob("*.tsv"):
        if "summary" in p.name.lower():
            return p

    print(
        f"[Richer_report] METABOLIC run completed but no summary TSV found in {out_dir}.",
        file=sys.stderr,
    )
    return None


def run_eggnog_mapper(
    annotation_dir: Path,
    out_dir: Path,
    threads: int,
) -> Optional[Path]:
    """
    Run eggNOG-mapper on all Prokka protein FASTAs (*.faa) and return
    the combined annotations TSV path if successful.
    """
    exe = find_executable(["emapper.py", "emapper"])
    if exe is None:
        print(
            "[Richer_report] eggNOG-mapper (emapper.py) not found on PATH. "
            "Skipping eggNOG run.",
            file=sys.stderr,
        )
        return None

    annotation_dir = annotation_dir.resolve()
    out_dir = out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    faa_files = sorted(annotation_dir.glob("*/*.faa"))
    if not faa_files:
        print(
            f"[Richer_report] No Prokka protein FASTAs (*.faa) found under {annotation_dir}; "
            "skipping eggNOG run.",
            file=sys.stderr,
        )
        return None

    combined_faa = out_dir / "eggnog_all_proteins.faa"
    with combined_faa.open("w") as out_fh:
        for faa in faa_files:
            with faa.open() as in_fh:
                out_fh.write(in_fh.read())

    prefix = out_dir / "eggnog"
    cmd = [
        exe,
        "-i",
        str(combined_faa),
        "-o",
        str(prefix),
        "--cpu",
        str(threads),
        "--override",
    ]

    data_dir = os.environ.get("EGGNOG_DATA_DIR")
    if data_dir:
        cmd.extend(["--data_dir", data_dir])

    db_path: Optional[Path] = None
    custom_db = os.environ.get("EGGNOG_DIAMOND_DB")
    if custom_db:
        db_path = Path(custom_db)
    elif data_dir:
        candidate = Path(data_dir) / "bacteria.dmnd"
        if candidate.exists():
            db_path = candidate
        else:
            default_db = Path(data_dir) / "eggnog_proteins.dmnd"
            if default_db.exists():
                db_path = default_db

    if db_path:
        cmd.extend(["--dmnd_db", str(db_path)])

    if not run_command(cmd):
        return None

    annotations = Path(f"{prefix}.emapper.annotations")
    if not annotations.exists():
        print(
            f"[Richer_report] eggNOG annotations file not found at {annotations}.",
            file=sys.stderr,
        )
        return None

    return annotations


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Enrich an existing CoSMIC LLM report (cosmic_llm_report.md) with "
            "additional functional context from DRAM, METABOLIC, HUMAnN, and "
            "eggNOG-mapper result tables."
        )
    )
    parser.add_argument(
        "--base-report",
        default="cosmic_llm_report.md",
        help="Path to the existing base report (default: cosmic_llm_report.md).",
    )
    parser.add_argument(
        "--mags-dir",
        default="Data",
        help="Directory containing MAG FASTA files (default: Data).",
    )
    parser.add_argument(
        "--annotation-dir",
        default="Annotation",
        help="Directory containing Prokka annotation subfolders (default: Annotation).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Number of CPU threads to use when running external tools (default: 8).",
    )
    parser.add_argument(
        "--run-dram",
        action="store_true",
        help="Run DRAM annotate+distill on MAGs before summarizing.",
    )
    parser.add_argument(
        "--run-metabolic",
        action="store_true",
        help="Run METABOLIC-G on MAGs before summarizing.",
    )
    parser.add_argument(
        "--run-eggnog",
        action="store_true",
        help="Run eggNOG-mapper on Prokka protein FASTAs before summarizing.",
    )
    parser.add_argument(
        "--dram-genome-summary",
        default=None,
        help="Optional DRAM genome summary TSV (per MAG). If not provided and "
        "--run-dram is set, this is auto-populated from DRAM output.",
    )
    parser.add_argument(
        "--dram-metabolism-summary",
        default=None,
        help="Optional DRAM metabolism summary TSV (per MAG). If not provided and "
        "--run-dram is set, this is auto-populated from DRAM output.",
    )
    parser.add_argument(
        "--metabolic-summary",
        default=None,
        help="Optional METABOLIC result summary TSV. If not provided and "
        "--run-metabolic is set, this is auto-populated from METABOLIC output.",
    )
    parser.add_argument(
        "--humann-pathway-abundance",
        default=None,
        help="Optional HUMAnN pathway abundance TSV.",
    )
    parser.add_argument(
        "--eggnog-annotations",
        default=None,
        help=(
            "Optional combined eggNOG-mapper annotations TSV. If not provided and "
            "--run-eggnog is set, this is auto-populated from eggNOG output."
        ),
    )
    parser.add_argument(
        "--output",
        default="cosmic_llm_rich_report.md",
        help="Output markdown file (default: cosmic_llm_rich_report.md).",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    base_path = Path(args.base_report)
    out_path = Path(args.output)

    mags_dir = Path(args.mags_dir)
    annotation_dir = Path(args.annotation_dir)

    # Optionally run external tools to generate functional summaries.
    if args.run_dram:
        dram_out = Path("DRAM_output")
        g_sum, m_sum = run_dram(mags_dir, dram_out, args.threads)
        if g_sum is not None:
            args.dram_genome_summary = str(g_sum)
        if m_sum is not None:
            args.dram_metabolism_summary = str(m_sum)

    if args.run_metabolic:
        metabolic_out = Path("METABOLIC_output")
        meta_sum = run_metabolic(mags_dir, metabolic_out, args.threads)
        if meta_sum is not None:
            args.metabolic_summary = str(meta_sum)

    if args.run_eggnog:
        eggnog_out = Path("eggNOG_output")
        eggnog_ann = run_eggnog_mapper(annotation_dir, eggnog_out, args.threads)
        if eggnog_ann is not None:
            args.eggnog_annotations = str(eggnog_ann)

    lines = read_text_lines(base_path)

    # Append additional sections if inputs are provided and exist
    if args.dram_genome_summary:
        lines.append("")
        lines.extend(summarize_dram_genome_summary(Path(args.dram_genome_summary)))

    if args.dram_metabolism_summary:
        lines.append("")
        lines.extend(summarize_dram_metabolism(Path(args.dram_metabolism_summary)))

    if args.metabolic_summary:
        lines.append("")
        lines.extend(summarize_metabolic_summary(Path(args.metabolic_summary)))

    if args.humann_pathway_abundance:
        lines.append("")
        lines.extend(summarize_humann_pathways(Path(args.humann_pathway_abundance)))

    if args.eggnog_annotations:
        lines.append("")
        lines.extend(summarize_eggnog_annotations(Path(args.eggnog_annotations)))

    with out_path.open("w") as fh:
        fh.write("\n".join(lines))

    print(f"Wrote enriched LLM report to {out_path}")


if __name__ == "__main__":
    main()
