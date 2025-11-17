"""
Optional script to enrich the CoSMIC LLM report with additional
functional context from external tools (DRAM, METABOLIC, HUMAnN,
eggNOG-mapper), if their result tables are available.

Usage (basic):

    source .venv/bin/activate
    python Richer_report.py \
        --base-report cosmic_llm_report.md \
        --output cosmic_llm_rich_report.md

You can also pass paths to additional tool outputs, e.g.:

    python Richer_report.py \
        --base-report cosmic_llm_report.md \
        --dram-genome-summary DRAM/genome_summary.tsv \
        --dram-metabolism-summary DRAM/metabolism_summary.tsv \
        --metabolic-summary METABOLIC/METABOLIC_result_summary.tsv \
        --humann-pathway-abundance humann/pathway_abundance.tsv \
        --eggnog-annotations eggnog/combined_annotations.tsv \
        --output cosmic_llm_rich_report.md

This script does NOT run those tools; it only integrates their
tabular outputs into a richer markdown report.
"""

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd


def read_text_lines(path: Path) -> List[str]:
    if not path.exists():
        raise SystemExit(f"Base report not found: {path}")
    with path.open() as fh:
        return fh.read().splitlines()


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

    df = pd.read_csv(path, sep="\t")
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

    df = pd.read_csv(path, sep="\t")

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
        "--dram-genome-summary",
        default=None,
        help="Optional DRAM genome summary TSV (per MAG).",
    )
    parser.add_argument(
        "--dram-metabolism-summary",
        default=None,
        help="Optional DRAM metabolism summary TSV (per MAG).",
    )
    parser.add_argument(
        "--metabolic-summary",
        default=None,
        help="Optional METABOLIC result summary TSV.",
    )
    parser.add_argument(
        "--humann-pathway-abundance",
        default=None,
        help="Optional HUMAnN pathway abundance TSV.",
    )
    parser.add_argument(
        "--eggnog-annotations",
        default=None,
        help="Optional combined eggNOG-mapper annotations TSV.",
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

