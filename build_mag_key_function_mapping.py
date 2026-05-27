#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


DEFAULT_RUNS = [
    ("Agelas", "V3V4", "run_Agelas_V3V4/cosmic_llm_report.md"),
    ("Agelas", "Full", "run_Agelas_full/cosmic_llm_report.md"),
    ("Oscarella", "V3V4", "run_Oscarella_V3V4/cosmic_llm_report.md"),
    ("Oscarella", "Full", "run_Oscarella_full/cosmic_llm_report.md"),
]


def parse_mag_annotation_section(report_path: Path):
    records = {}
    in_species_section = False
    current_mag = None
    current_list = None

    with report_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()

            if stripped == "## Species / MAG Annotations and Relative Abundances":
                in_species_section = True
                continue

            if in_species_section and stripped.startswith("## ") and stripped != "## Species / MAG Annotations and Relative Abundances":
                break

            if not in_species_section:
                continue

            if stripped.startswith("### MAG "):
                current_mag = stripped.replace("### MAG ", "", 1)
                records[current_mag] = {
                    "annotated_cds_count": "",
                    "non_hypothetical_cds_count": "",
                    "products": [],
                    "ecs": [],
                    "cogs": [],
                }
                current_list = None
                continue

            if current_mag is None:
                continue

            if stripped.startswith("- Annotated CDS count:"):
                records[current_mag]["annotated_cds_count"] = stripped.split(":", 1)[1].strip()
                current_list = None
                continue

            if stripped.startswith("- Non-hypothetical CDS count:"):
                records[current_mag]["non_hypothetical_cds_count"] = stripped.split(":", 1)[1].strip()
                current_list = None
                continue

            if stripped == "- Most frequent annotated products:":
                current_list = "products"
                continue

            if stripped == "- Top EC numbers (by CDS count):":
                current_list = "ecs"
                continue

            if stripped == "- Top COGs (by CDS count):":
                current_list = "cogs"
                continue

            if stripped.startswith("### MAG ") or stripped.startswith("## "):
                current_list = None
                continue

            if stripped.startswith("- ") and current_list:
                item = stripped[2:].strip()
                records[current_mag][current_list].append(item)
                continue

            if stripped == "":
                current_list = current_list

    return records


def classify_presence(v3v4_mags, full_mags):
    classes = {}
    all_mags = v3v4_mags | full_mags
    for mag in all_mags:
        if mag in v3v4_mags and mag in full_mags:
            classes[mag] = "both"
        elif mag in v3v4_mags:
            classes[mag] = "v3v4_only"
        else:
            classes[mag] = "full_only"
    return classes


def top_items(values, top_n):
    return " | ".join(values[:top_n])


def write_tsv(path: Path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Build MAG key-function mapping tables from CoSMIC Sidekick report files."
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="How many top products/ECs/COGs to keep in key-function columns (default: 5).",
    )
    parser.add_argument(
        "--output-dir",
        default=".",
        help="Directory for output TSV files (default: current directory).",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    parsed = {}
    for dataset, run_type, report in DEFAULT_RUNS:
        report_path = Path(report).resolve()
        if not report_path.exists():
            raise FileNotFoundError(f"Report not found: {report_path}")
        parsed.setdefault(dataset, {})[run_type] = {
            "report_path": str(report_path),
            "records": parse_mag_annotation_section(report_path),
        }

    long_rows = []
    comparison_rows = []

    for dataset in sorted(parsed):
        v3v4 = parsed[dataset].get("V3V4", {"records": {}})["records"]
        full = parsed[dataset].get("Full", {"records": {}})["records"]
        classes = classify_presence(set(v3v4.keys()), set(full.keys()))

        for run_type in ("V3V4", "Full"):
            run_data = parsed[dataset].get(run_type)
            if not run_data:
                continue
            report_path = run_data["report_path"]
            for mag in sorted(run_data["records"]):
                rec = run_data["records"][mag]
                long_rows.append(
                    {
                        "dataset": dataset,
                        "run_type": run_type,
                        "mag_id": mag,
                        "presence_class": classes.get(mag, ""),
                        "annotated_cds_count": rec["annotated_cds_count"],
                        "non_hypothetical_cds_count": rec["non_hypothetical_cds_count"],
                        f"key_products_top{args.top_n}": top_items(rec["products"], args.top_n),
                        f"key_ec_top{args.top_n}": top_items(rec["ecs"], args.top_n),
                        f"key_cog_top{args.top_n}": top_items(rec["cogs"], args.top_n),
                        "source_report": report_path,
                    }
                )

        for mag in sorted(set(v3v4.keys()) | set(full.keys())):
            v = v3v4.get(
                mag,
                {"annotated_cds_count": "", "non_hypothetical_cds_count": "", "products": [], "ecs": [], "cogs": []},
            )
            f = full.get(
                mag,
                {"annotated_cds_count": "", "non_hypothetical_cds_count": "", "products": [], "ecs": [], "cogs": []},
            )
            comparison_rows.append(
                {
                    "dataset": dataset,
                    "mag_id": mag,
                    "presence_class": classes.get(mag, ""),
                    "v3v4_annotated_cds_count": v["annotated_cds_count"],
                    "v3v4_non_hypothetical_cds_count": v["non_hypothetical_cds_count"],
                    f"v3v4_key_products_top{args.top_n}": top_items(v["products"], args.top_n),
                    f"v3v4_key_ec_top{args.top_n}": top_items(v["ecs"], args.top_n),
                    f"v3v4_key_cog_top{args.top_n}": top_items(v["cogs"], args.top_n),
                    "full_annotated_cds_count": f["annotated_cds_count"],
                    "full_non_hypothetical_cds_count": f["non_hypothetical_cds_count"],
                    f"full_key_products_top{args.top_n}": top_items(f["products"], args.top_n),
                    f"full_key_ec_top{args.top_n}": top_items(f["ecs"], args.top_n),
                    f"full_key_cog_top{args.top_n}": top_items(f["cogs"], args.top_n),
                }
            )

    long_path = output_dir / "Agelas_Oscarella_MAG_key_function_mapping.tsv"
    comparison_path = output_dir / "Agelas_Oscarella_MAG_key_function_comparison.tsv"

    write_tsv(
        long_path,
        [
            "dataset",
            "run_type",
            "mag_id",
            "presence_class",
            "annotated_cds_count",
            "non_hypothetical_cds_count",
            f"key_products_top{args.top_n}",
            f"key_ec_top{args.top_n}",
            f"key_cog_top{args.top_n}",
            "source_report",
        ],
        long_rows,
    )

    write_tsv(
        comparison_path,
        [
            "dataset",
            "mag_id",
            "presence_class",
            "v3v4_annotated_cds_count",
            "v3v4_non_hypothetical_cds_count",
            f"v3v4_key_products_top{args.top_n}",
            f"v3v4_key_ec_top{args.top_n}",
            f"v3v4_key_cog_top{args.top_n}",
            "full_annotated_cds_count",
            "full_non_hypothetical_cds_count",
            f"full_key_products_top{args.top_n}",
            f"full_key_ec_top{args.top_n}",
            f"full_key_cog_top{args.top_n}",
        ],
        comparison_rows,
    )

    print(f"Wrote: {long_path}")
    print(f"Wrote: {comparison_path}")
    print(f"Rows (long): {len(long_rows)}")
    print(f"Rows (comparison): {len(comparison_rows)}")


if __name__ == "__main__":
    main()
