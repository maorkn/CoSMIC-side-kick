import argparse
import gzip
import os
import random
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd
import yaml
from Bio import SeqIO, pairwise2


@dataclass
class RRnaRecord:
    rrna_uid: str
    mag_id: str
    mag_fasta_path: str
    sequence_header: str
    rrna_id: str
    rrna_type: str
    start: int
    end: int
    strand: str
    length: int
    barrnap_kingdom: str
    product: str
    sequence: str


def load_config(path: Optional[str]) -> Dict:
    cfg_path = Path(path) if path else Path("config.yaml")
    if not cfg_path.exists():
        raise SystemExit(f"Config file not found: {cfg_path}")
    with cfg_path.open() as fh:
        return yaml.safe_load(fh) or {}


def get_mag_id(path: Path) -> str:
    name = path.name
    for ext in [".fasta.gz", ".fa.gz", ".fna.gz", ".fasta", ".fa", ".fna"]:
        if name.endswith(ext):
            return name[: -len(ext)]
    return path.stem


def open_fasta_for_seqio(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def run_barrnap_on_mag(
    fasta_path: Path, kingdoms: Sequence[str]
) -> List[Tuple[str, str, int, int, str, str, str, str]]:
    """
    Run barrnap on a MAG fasta for the given kingdoms.

    Returns list of tuples:
    (kingdom, seqid, rrna_id, rrna_type, start, end, strand, product)
    with start/end as ints.
    """
    from shutil import which

    if which("barrnap") is None:
        raise SystemExit(
            "barrnap not found on PATH. Please install barrnap "
            "and ensure it is available in your environment."
        )

    results: List[Tuple[str, str, str, str, int, int, str, str]] = []

    # Decompress to a temporary uncompressed fasta for barrnap
    import tempfile

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as tmp_fa:
        with open_fasta_for_seqio(fasta_path) as handle:
            tmp_fa.write(handle.read())
        tmp_path = Path(tmp_fa.name)

    try:
        for kingdom in kingdoms:
            cmd = ["barrnap", "--kingdom", kingdom, str(tmp_path)]
            proc = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
            gff_text = proc.stdout.splitlines()
            for line in gff_text:
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) != 9:
                    continue
                seqid, source, feature_type, start, end, score, strand, phase, attributes = (
                    parts
                )
                if feature_type != "rRNA":
                    continue

                attrs: Dict[str, str] = {}
                for field in attributes.split(";"):
                    if not field:
                        continue
                    if "=" in field:
                        key, val = field.split("=", 1)
                        attrs[key] = val

                product = attrs.get("product", attrs.get("Name", ""))
                rrna_id = attrs.get("ID", "")

                rrna_type = ""
                product_upper = product.upper()
                if "16S" in product_upper:
                    rrna_type = "16S"
                elif "18S" in product_upper:
                    rrna_type = "18S"
                else:
                    # Skip other rRNAs (5S, 23S, etc.)
                    continue

                results.append(
                    (
                        kingdom,
                        seqid,
                        rrna_id,
                        rrna_type,
                        int(start),
                        int(end),
                        strand,
                        product,
                    )
                )
    finally:
        try:
            tmp_path.unlink()
        except Exception:
            pass

    return results


def extract_rrna_from_mags(
    mags_dir: Path, kingdoms: Sequence[str]
) -> List[RRnaRecord]:
    mags: List[Path] = sorted(
        p for p in mags_dir.iterdir() if p.suffix in {".fa", ".fasta", ".fna", ".gz"}
    )
    if not mags:
        raise SystemExit(f"No MAG fasta files found in {mags_dir}")

    all_records: List[RRnaRecord] = []

    for mag_path in mags:
        mag_id = get_mag_id(mag_path)

        # Load sequences into memory (dict of id -> sequence)
        seq_records: Dict[str, str] = {}
        with open_fasta_for_seqio(mag_path) as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                seq_records[rec.id] = str(rec.seq).upper()

        gff_hits = run_barrnap_on_mag(mag_path, kingdoms)
        for kingdom, seqid, rrna_id, rrna_type, start, end, strand, product in gff_hits:
            # Barrnap coordinates are 1-based inclusive
            seq = seq_records.get(seqid)
            if not seq:
                continue
            subseq = seq[start - 1 : end]
            if strand == "-":
                # Reverse complement
                complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
                subseq = subseq.translate(complement)[::-1]

            length = len(subseq)
            rrna_uid = f"{mag_id}|{seqid}|{rrna_type}|{start}-{end}"

            record = RRnaRecord(
                rrna_uid=rrna_uid,
                mag_id=mag_id,
                mag_fasta_path=str(mag_path),
                sequence_header=seqid,
                rrna_id=rrna_id,
                rrna_type=rrna_type,
                start=start,
                end=end,
                strand=strand,
                length=length,
                barrnap_kingdom=kingdom,
                product=product,
                sequence=subseq,
            )
            all_records.append(record)

    return all_records


def save_rrna_outputs(records: List[RRnaRecord], out_prefix: str = "barrnap_rrna"):
    if not records:
        print("No 16S/18S rRNA features found in MAGs.", file=sys.stderr)
        return

    df = pd.DataFrame(
        [
            {
                "rrna_uid": r.rrna_uid,
                "mag_id": r.mag_id,
                "mag_fasta_path": r.mag_fasta_path,
                "sequence_header": r.sequence_header,
                "rrna_id": r.rrna_id,
                "rrna_type": r.rrna_type,
                "start": r.start,
                "end": r.end,
                "strand": r.strand,
                "length": r.length,
                "barrnap_kingdom": r.barrnap_kingdom,
                "product": r.product,
                "sequence": r.sequence,
            }
            for r in records
        ]
    )
    csv_path = Path(f"{out_prefix}_mapping.csv")
    df.to_csv(csv_path, index=False)

    fasta_path = Path(f"{out_prefix}_sequences.fasta")
    with fasta_path.open("w") as fh:
        for r in records:
            fh.write(f">{r.rrna_uid}\n")
            fh.write(f"{r.sequence}\n")

    print(f"Saved rRNA mapping CSV to {csv_path}")
    print(f"Saved rRNA FASTA to {fasta_path}")


def load_metabarcoding_table(
    path: Path,
    id_column: str,
    sequence_column: str,
    abundance_columns: Optional[List[str]] = None,
) -> Tuple[pd.DataFrame, List[str]]:
    df = pd.read_csv(path)

    if id_column not in df.columns:
        raise SystemExit(f"ID column '{id_column}' not found in {path}")
    if sequence_column not in df.columns:
        raise SystemExit(f"Sequence column '{sequence_column}' not found in {path}")

    if abundance_columns is None:
        # Auto-detect numeric columns except id and sequence
        candidates = [
            c
            for c in df.columns
            if c not in {id_column, sequence_column} and pd.api.types.is_numeric_dtype(df[c])
        ]
        abundance_columns = candidates
    else:
        for col in abundance_columns:
            if col not in df.columns:
                raise SystemExit(
                    f"Configured abundance column '{col}' not found in {path}"
                )

    return df, abundance_columns


def compute_identity(seq1: str, seq2: str) -> float:
    s1 = seq1.upper().replace("U", "T")
    s2 = seq2.upper().replace("U", "T")
    if not s1 or not s2:
        return 0.0

    # Fast path: equal length, no gaps
    if len(s1) == len(s2):
        matches = sum(1 for a, b in zip(s1, s2) if a == b and a != "N" and b != "N")
        length = sum(1 for a, b in zip(s1, s2) if a != "N" and b != "N")
        if length == 0:
            return 0.0
        return matches / length

    # Fallback: global alignment
    # Skip clearly incompatible lengths
    max_len = max(len(s1), len(s2))
    if abs(len(s1) - len(s2)) > max(100, int(0.2 * max_len)):
        return 0.0

    alignment = pairwise2.align.globalxx(s1, s2, one_alignment_only=True)[0]
    aligned1, aligned2, score, start, end = alignment
    length = len(aligned1)
    if length == 0:
        return 0.0
    matches = sum(
        1 for a, b in zip(aligned1, aligned2) if a == b and a != "-" and b != "-"
    )
    return matches / length


def map_metabarcoding_to_rrna(
    rrna_records: List[RRnaRecord],
    meta_df: pd.DataFrame,
    id_col: str,
    seq_col: str,
    abundance_cols: List[str],
    identity_threshold: float,
) -> pd.DataFrame:
    rrna_by_uid = {r.rrna_uid: r for r in rrna_records}
    rrna_items = list(rrna_by_uid.items())

    mappings: List[Dict] = []

    for _, row in meta_df.iterrows():
        meta_id = row[id_col]
        meta_seq = str(row[seq_col])
        if not isinstance(meta_seq, str) or not meta_seq.strip():
            continue
        meta_seq = meta_seq.strip()

        for rrna_uid, r in rrna_items:
            pid = compute_identity(meta_seq, r.sequence)
            if pid >= identity_threshold:
                entry: Dict = {
                    "metabarcoding_id": meta_id,
                    "metabarcoding_sequence": meta_seq,
                    "rrna_uid": rrna_uid,
                    "mag_id": r.mag_id,
                    "sequence_header": r.sequence_header,
                    "rrna_type": r.rrna_type,
                    "identity": pid,
                }
                for col in abundance_cols:
                    entry[col] = row[col]
                mappings.append(entry)

    if not mappings:
        print(
            "No metabarcoding sequences mapped to MAG rRNA at the given threshold.",
            file=sys.stderr,
        )
        return pd.DataFrame()

    return pd.DataFrame(mappings)


def run_prokka_on_mags(
    mags_by_id: Dict[str, Path],
    mag_ids_to_annotate: Sequence[str],
    out_dir: Path,
):
    from shutil import which

    if which("prokka") is None:
        print(
            "Warning: prokka not found on PATH. Skipping annotation.",
            file=sys.stderr,
        )
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    import tempfile

    for mag_id in sorted(set(mag_ids_to_annotate)):
        mag_path = mags_by_id.get(mag_id)
        if mag_path is None:
            continue
        mag_out = out_dir / mag_id
        if mag_out.exists():
            # Assume already annotated
            print(f"Annotation already exists for {mag_id} at {mag_out}, skipping.")
            continue

        # Prokka prefers uncompressed fasta
        if str(mag_path).endswith(".gz"):
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fa", delete=False
            ) as tmp_fa:
                with open_fasta_for_seqio(mag_path) as handle:
                    tmp_fa.write(handle.read())
                tmp_in = Path(tmp_fa.name)
            input_path = tmp_in
        else:
            input_path = mag_path
            tmp_in = None

        cmd = [
            "prokka",
            "--outdir",
            str(mag_out),
            "--prefix",
            mag_id,
            str(input_path),
        ]
        print(f"Running: {' '.join(cmd)}")
        env = os.environ.copy()
        if "PROKKA_DBDIR" not in env:
            default_dbdir = "/usr/share/prokka/db"
            if Path(default_dbdir).exists():
                env["PROKKA_DBDIR"] = default_dbdir
        try:
            subprocess.run(cmd, check=True, env=env)
        finally:
            if tmp_in is not None:
                try:
                    tmp_in.unlink()
                except Exception:
                    pass


def run_pipeline(args: argparse.Namespace):
    config = load_config(args.config)

    mags_dir = Path(args.mags_dir or config.get("mags_dir", "Data"))
    metabarcoding_csv = Path(
        args.metabarcoding or config.get("metabarcoding_csv", "metabarcoding.csv")
    )
    id_col = config.get("metabarcoding_id_column", "id")
    seq_col = config.get("metabarcoding_sequence_column", "sequence")
    abundance_cols = config.get("metabarcoding_abundance_columns")
    identity_threshold = float(config.get("identity_threshold", 0.97))
    kingdoms = config.get("barrnap_kingdoms", ["bac", "arc", "euk"])
    annotation_tool = config.get("annotation_tool", "prokka")
    annotation_output_dir = Path(config.get("annotation_output_dir", "Annotation"))

    print(f"Using MAGs directory: {mags_dir}")
    mags_dir = mags_dir.resolve()

    rrna_records = extract_rrna_from_mags(mags_dir, kingdoms)
    save_rrna_outputs(rrna_records)

    if not metabarcoding_csv.exists():
        print(
            f"Metabarcoding CSV '{metabarcoding_csv}' not found. "
            "Pipeline will stop after rRNA extraction.",
            file=sys.stderr,
        )
        return

    meta_df, abundance_cols = load_metabarcoding_table(
        metabarcoding_csv, id_col, seq_col, abundance_cols
    )

    mapping_df = map_metabarcoding_to_rrna(
        rrna_records,
        meta_df,
        id_col,
        seq_col,
        abundance_cols,
        identity_threshold,
    )

    if mapping_df.empty:
        print("No mappings found; skipping annotation.", file=sys.stderr)
        return

    mapping_csv = Path("metabarcoding_to_MAG_mapping.csv")
    mapping_df.to_csv(mapping_csv, index=False)
    print(f"Saved metabarcoding-to-MAG mapping to {mapping_csv}")

    # Annotation
    if annotation_tool.lower() == "prokka":
        mags_by_id: Dict[str, Path] = {}
        for p in mags_dir.iterdir():
            if p.suffix in {".fa", ".fasta", ".fna", ".gz"}:
                mags_by_id[get_mag_id(p)] = p
        mag_ids_to_annotate = mapping_df["mag_id"].unique().tolist()
        run_prokka_on_mags(mags_by_id, mag_ids_to_annotate, annotation_output_dir)
    else:
        print(
            f"Annotation tool '{annotation_tool}' is not implemented in this script.",
            file=sys.stderr,
        )


def cmd_extract_rrna(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    mags_dir = Path(args.mags_dir or config.get("mags_dir", "Data"))
    kingdoms = config.get("barrnap_kingdoms", ["bac", "arc", "euk"])

    print(f"Using MAGs directory: {mags_dir}")
    mags_dir = mags_dir.resolve()

    rrna_records = extract_rrna_from_mags(mags_dir, kingdoms)
    save_rrna_outputs(rrna_records)


def cmd_annotate_test(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    mags_dir = Path(args.mags_dir or config.get("mags_dir", "Data"))
    mags_dir = mags_dir.resolve()
    annotation_output_dir = Path(config.get("annotation_output_dir", "Annotation"))

    mapping_csv = Path(args.rrna_mapping or "barrnap_rrna_mapping.csv")
    if not mapping_csv.exists():
        raise SystemExit(
            f"rRNA mapping CSV not found: {mapping_csv}. "
            "Run 'extract-rrna' first to generate it."
        )

    df = pd.read_csv(mapping_csv)
    if "mag_id" not in df.columns:
        raise SystemExit(
            f"Column 'mag_id' not found in {mapping_csv}; "
            "cannot determine which MAGs have rRNA hits."
        )

    mag_ids = sorted(df["mag_id"].dropna().unique().tolist())
    if not mag_ids:
        print(
            "No MAG IDs found in rRNA mapping CSV; nothing to annotate.",
            file=sys.stderr,
        )
        return

    num = max(1, int(args.num_mags))
    num = min(num, len(mag_ids))
    chosen_mag_ids = random.sample(mag_ids, num)

    print(f"Available MAGs with rRNA hits: {mag_ids}")
    print(f"Randomly selected MAGs for Prokka annotation (n={num}): {chosen_mag_ids}")

    mags_by_id: Dict[str, Path] = {}
    for p in mags_dir.iterdir():
        if p.suffix in {".fa", ".fasta", ".fna", ".gz"}:
            mags_by_id[get_mag_id(p)] = p

    run_prokka_on_mags(mags_by_id, chosen_mag_ids, annotation_output_dir)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Pipeline to link MAGs with metabarcoding 16S/18S and annotate MAGs."
    )
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser(
        "run",
        help="Run the full pipeline (rRNA extraction, mapping, annotation).",
    )
    run_parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file (default: config.yaml).",
    )
    run_parser.add_argument(
        "--mags-dir",
        default=None,
        help="Directory containing MAG fasta(.gz) files (overrides config).",
    )
    run_parser.add_argument(
        "--metabarcoding",
        default=None,
        help="Metabarcoding CSV file (overrides config).",
    )

    extract_parser = subparsers.add_parser(
        "extract-rrna",
        help="Run only rRNA extraction with Barrnap and write mapping/FASTA files.",
    )
    extract_parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file (default: config.yaml).",
    )
    extract_parser.add_argument(
        "--mags-dir",
        default=None,
        help="Directory containing MAG fasta(.gz) files (overrides config).",
    )

    annotate_test_parser = subparsers.add_parser(
        "annotate-test",
        help=(
            "Run Prokka on a random subset of MAGs that have rRNA hits "
            "according to an existing mapping file."
        ),
    )
    annotate_test_parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file (default: config.yaml).",
    )
    annotate_test_parser.add_argument(
        "--mags-dir",
        default=None,
        help="Directory containing MAG fasta(.gz) files (overrides config).",
    )
    annotate_test_parser.add_argument(
        "--rrna-mapping",
        default="barrnap_rrna_mapping.csv",
        help="Path to existing rRNA mapping CSV (default: barrnap_rrna_mapping.csv).",
    )
    annotate_test_parser.add_argument(
        "--num-mags",
        type=int,
        default=5,
        help="How many MAGs to annotate (default: 5).",
    )

    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "run":
        run_pipeline(args)
    elif args.command == "extract-rrna":
        cmd_extract_rrna(args)
    elif args.command == "annotate-test":
        cmd_annotate_test(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
