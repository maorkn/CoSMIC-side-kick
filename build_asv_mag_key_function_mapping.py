#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path
from urllib.parse import unquote


def clean_text(text: str) -> str:
    cleaned = (
        text.replace("\\_", "_")
        .replace("\\[", "[")
        .replace("\\]", "]")
        .strip()
    )
    return unquote(cleaned)


def parse_asv_report(report_path: Path):
    asv_records = []

    current_asv_id = None
    current_seq_len_bp = ""
    current_abund = {}
    current_hit = None
    current_list = None

    def finalize_hit():
        nonlocal current_hit
        if current_hit is None:
            return
        asv_records.append(current_hit)
        current_hit = None

    with report_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()

            m_asv = re.match(r"^### Metabarcoding ID\s+(.+)$", stripped)
            if m_asv:
                finalize_hit()
                current_asv_id = clean_text(m_asv.group(1))
                current_seq_len_bp = ""
                current_abund = {}
                current_list = None
                continue

            if current_asv_id is None:
                continue

            m_seq = re.match(r"^\*\s+Sequence length:\s+(\d+)\s+bp$", stripped)
            if m_seq:
                current_seq_len_bp = m_seq.group(1)
                continue

            m_ab = re.match(r"^\*\s+(.+?):\s+([0-9.]+)$", stripped)
            if m_ab and not stripped.startswith("* Sequence length:") and "identity" not in stripped:
                key = clean_text(m_ab.group(1))
                val = m_ab.group(2)
                current_abund[key] = val
                continue

            m_hit = re.match(
                r"^\*\s+MAG\s+(.+?),\s+contig\s+(.+?),\s+rRNA hit\s+(.+?),\s+identity\s+([0-9.]+)%$",
                stripped,
            )
            if m_hit:
                finalize_hit()
                current_hit = {
                    "metabarcoding_id": current_asv_id,
                    "sequence_length_bp": current_seq_len_bp,
                    "mag_id": clean_text(m_hit.group(1)),
                    "contig": clean_text(m_hit.group(2)),
                    "rrna_hit": clean_text(m_hit.group(3)),
                    "identity_percent": m_hit.group(4),
                    "function_source": "",
                    "key_functions_top5": "",
                    "mag_products_top5": "",
                    "mag_ec_top5": "",
                    "mag_cog_top5": "",
                    "heatwave - dying": current_abund.get("heatwave - dying", ""),
                    "heatwave-no sings of necrosis": current_abund.get("heatwave-no sings of necrosis", ""),
                    "pre-heatwave- healthy": current_abund.get("pre-heatwave- healthy", ""),
                    "post heatwave- healthy": current_abund.get("post heatwave- healthy", ""),
                    "heatwave-showing sings of necrosis": current_abund.get("heatwave-showing sings of necrosis", ""),
                    "_contig_highlights": [],
                    "_mag_products": [],
                    "_mag_ecs": [],
                    "_mag_cogs": [],
                }
                current_list = None
                continue

            if current_hit is None:
                continue

            if stripped == "* Contig annotation highlights:":
                current_list = "_contig_highlights"
                continue
            if stripped == "* MAG products (top):":
                current_list = "_mag_products"
                continue
            if stripped == "* MAG EC numbers (top):":
                current_list = "_mag_ecs"
                continue
            if stripped == "* MAG COGs (top):":
                current_list = "_mag_cogs"
                continue

            if current_list and stripped.startswith("* "):
                current_hit[current_list].append(clean_text(stripped[2:].strip()))
                continue

            if stripped.startswith("### Metabarcoding ID "):
                current_list = None

    finalize_hit()

    for rec in asv_records:
        contig_top = rec["_contig_highlights"][:5]
        mag_prod_top = rec["_mag_products"][:5]
        if contig_top:
            rec["function_source"] = "contig_annotation_highlights"
            rec["key_functions_top5"] = " | ".join(contig_top)
        else:
            rec["function_source"] = "mag_products_top"
            rec["key_functions_top5"] = " | ".join(mag_prod_top)
        rec["mag_products_top5"] = " | ".join(mag_prod_top)
        rec["mag_ec_top5"] = " | ".join(rec["_mag_ecs"][:5])
        rec["mag_cog_top5"] = " | ".join(rec["_mag_cogs"][:5])
        for k in ["_contig_highlights", "_mag_products", "_mag_ecs", "_mag_cogs"]:
            rec.pop(k, None)

    return asv_records


def write_tsv(path: Path, rows):
    fields = [
        "metabarcoding_id",
        "sequence_length_bp",
        "mag_id",
        "contig",
        "rrna_hit",
        "identity_percent",
        "function_source",
        "key_functions_top5",
        "mag_products_top5",
        "mag_ec_top5",
        "mag_cog_top5",
        "heatwave - dying",
        "heatwave-no sings of necrosis",
        "pre-heatwave- healthy",
        "post heatwave- healthy",
        "heatwave-showing sings of necrosis",
        "source_report",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build ASV->MAG key-function mapping from ASV-centric CoSMIC markdown report."
    )
    parser.add_argument("--report", required=True, help="Path to cosmic_llm_report_*.md ASV report.")
    parser.add_argument("--output", required=True, help="Output TSV path.")
    args = parser.parse_args()

    report = Path(args.report).resolve()
    out = Path(args.output).resolve()

    rows = parse_asv_report(report)
    for r in rows:
        r["source_report"] = str(report)

    write_tsv(out, rows)
    print(f"Wrote: {out}")
    print(f"Rows: {len(rows)}")


if __name__ == "__main__":
    main()
