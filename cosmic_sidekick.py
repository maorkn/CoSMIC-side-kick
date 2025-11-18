import argparse
import gzip
import os
import random
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Set

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
        print(f"[CoSMIC] Running Barrnap on MAG {mag_id} ({mag_path.name})...")

        # Load sequences into memory (dict of id -> sequence)
        seq_records: Dict[str, str] = {}
        with open_fasta_for_seqio(mag_path) as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                seq_records[rec.id] = str(rec.seq).upper()

        gff_hits = run_barrnap_on_mag(mag_path, kingdoms)
        if not gff_hits:
            print(f"[CoSMIC]   Barrnap returned no 16S/18S hits for {mag_id}.")
            continue

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

        n_16s = sum(1 for _, _, _, rrna_type, *_ in gff_hits if rrna_type == "16S")
        n_18s = sum(1 for _, _, _, rrna_type, *_ in gff_hits if rrna_type == "18S")
        print(
            f"[CoSMIC]   Barrnap hits for {mag_id}: {len(gff_hits)} total "
            f"({n_16s} × 16S, {n_18s} × 18S)"
        )

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
    mag_to_contigs: Optional[Dict[str, Set[str]]] = None,
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

        # If a contig filter is provided for this MAG, restrict annotation
        # to only those contigs (saves time and focuses on CoSMIC-linked regions).
        contig_filter: Optional[Set[str]] = None
        if mag_to_contigs is not None:
            contig_filter = mag_to_contigs.get(mag_id)

        tmp_in: Optional[Path] = None

        if contig_filter:
            # Build temporary FASTA containing only the selected contigs.
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fa", delete=False
            ) as tmp_fa:
                tmp_in = Path(tmp_fa.name)
                kept_any = False
                with open_fasta_for_seqio(mag_path) as handle:
                    for rec in SeqIO.parse(handle, "fasta"):
                        if rec.id in contig_filter:
                            SeqIO.write(rec, tmp_fa, "fasta")
                            kept_any = True
            if not kept_any:
                # No contigs to annotate for this MAG under this filter.
                if tmp_in is not None:
                    try:
                        tmp_in.unlink()
                    except Exception:
                        pass
                print(
                    f"No contigs with CoSMIC-linked rRNA hits found for MAG {mag_id}; skipping.",
                    file=sys.stderr,
                )
                continue
            input_path = tmp_in
        else:
            # Prokka prefers uncompressed fasta; if needed, write a temporary copy.
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
    output_dir = Path(args.output_dir or config.get("output_dir", "."))
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    metabarcoding_csv = Path(
        args.metabarcoding or config.get("metabarcoding_csv", "metabarcoding.csv")
    )
    metabarcoding_csv = metabarcoding_csv.resolve()
    id_col = config.get("metabarcoding_id_column", "id")
    seq_col = config.get("metabarcoding_sequence_column", "sequence")
    abundance_cols = config.get("metabarcoding_abundance_columns")
    identity_threshold = float(config.get("identity_threshold", 0.97))
    kingdoms = config.get("barrnap_kingdoms", ["bac", "arc", "euk"])
    annotation_tool = config.get("annotation_tool", "prokka")
    annotation_output_cfg = Path(config.get("annotation_output_dir", "Annotation"))
    if annotation_output_cfg.is_absolute():
        annotation_output_dir = annotation_output_cfg
    else:
        annotation_output_dir = (output_dir / annotation_output_cfg).resolve()
    annotation_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[CoSMIC] Using MAGs directory: {mags_dir}")
    print(f"[CoSMIC] Output directory: {output_dir}")
    print(f"[CoSMIC] Annotation output directory: {annotation_output_dir}")
    print(f"[CoSMIC] Metabarcoding table: {metabarcoding_csv}")
    print(f"[CoSMIC] Identity threshold: {identity_threshold:.2%}")
    mags_dir = mags_dir.resolve()

    rrna_records = extract_rrna_from_mags(mags_dir, kingdoms)
    save_rrna_outputs(rrna_records, out_prefix=str(output_dir / "barrnap_rrna"))

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
    print(
        f"[CoSMIC] Loaded metabarcoding table with {len(meta_df)} rows "
        f"and {len(abundance_cols)} abundance columns."
    )

    mapping_df = map_metabarcoding_to_rrna(
        rrna_records,
        meta_df,
        id_col,
        seq_col,
        abundance_cols,
        identity_threshold,
    )
    print(
        f"[CoSMIC] Mapping complete: {len(mapping_df)} CoSMIC hits "
        f"covering {mapping_df['mag_id'].nunique() if not mapping_df.empty else 0} MAGs."
    )

    if mapping_df.empty:
        print("No mappings found; skipping annotation.", file=sys.stderr)
        return

    mapping_csv = output_dir / "metabarcoding_to_MAG_mapping.csv"
    mapping_df.to_csv(mapping_csv, index=False)
    print(f"Saved metabarcoding-to-MAG mapping to {mapping_csv}")

    # Annotation
    if annotation_tool.lower() == "prokka":
        mags_by_id: Dict[str, Path] = {}
        for p in mags_dir.iterdir():
            if p.suffix in {".fa", ".fasta", ".fna", ".gz"}:
                mags_by_id[get_mag_id(p)] = p
        mag_ids_to_annotate = mapping_df["mag_id"].unique().tolist()
        print(f"[CoSMIC] Annotating {len(mag_ids_to_annotate)} MAG(s) with Prokka...")

        # Determine which contigs within each MAG have CoSMIC-linked rRNA hits.
        mag_to_contigs: Dict[str, Set[str]] = {}
        if "sequence_header" in mapping_df.columns:
            for mag_id, group in mapping_df.groupby("mag_id"):
                contigs = set(group["sequence_header"].dropna().astype(str))
                if contigs:
                    mag_to_contigs[str(mag_id)] = contigs

        run_prokka_on_mags(
            mags_by_id,
            mag_ids_to_annotate,
            annotation_output_dir,
            mag_to_contigs=mag_to_contigs if mag_to_contigs else None,
        )
    else:
        print(
            f"Annotation tool '{annotation_tool}' is not implemented in this script.",
            file=sys.stderr,
        )


def cmd_extract_rrna(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    mags_dir = Path(args.mags_dir or config.get("mags_dir", "Data"))
    kingdoms = config.get("barrnap_kingdoms", ["bac", "arc", "euk"])
    output_dir = Path(args.output_dir or config.get("output_dir", "."))

    print(f"Using MAGs directory: {mags_dir}")
    mags_dir = mags_dir.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    rrna_records = extract_rrna_from_mags(mags_dir, kingdoms)

    # Temporarily change working directory for saving outputs
    cwd = Path.cwd()
    try:
        os.chdir(output_dir)
        save_rrna_outputs(rrna_records)
    finally:
        os.chdir(cwd)


def cmd_annotate_test(args: argparse.Namespace) -> None:
    config = load_config(args.config)
    mags_dir = Path(args.mags_dir or config.get("mags_dir", "Data"))
    mags_dir = mags_dir.resolve()

    output_dir = Path(args.output_dir or config.get("output_dir", ".")).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    annotation_output_dir = output_dir / config.get("annotation_output_dir", "Annotation")

    if args.rrna_mapping:
        mapping_csv = Path(args.rrna_mapping)
    else:
        mapping_csv = output_dir / "barrnap_rrna_mapping.csv"
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


def cmd_report(args: argparse.Namespace) -> None:
    """
    Generate a markdown-style report summarizing MAG metadata,
    CoSMIC experiment metadata, and species-level annotations +
    relative abundances, and append an LLM query prompt.
    """
    config = load_config(args.config)

    output_dir = Path(args.output_dir or config.get("output_dir", ".")).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.rrna_mapping:
        rrna_mapping_csv = Path(args.rrna_mapping)
    else:
        rrna_mapping_csv = output_dir / "barrnap_rrna_mapping.csv"

    if args.annotation_dir:
        mag_annotation_dir = Path(args.annotation_dir)
    else:
        mag_annotation_dir = output_dir / config.get("annotation_output_dir", "Annotation")

    if args.metabarcoding_mapping:
        meta_mapping_csv = Path(args.metabarcoding_mapping)
    else:
        meta_mapping_csv = output_dir / "metabarcoding_to_MAG_mapping.csv"
    experiment_metadata_yaml = config.get("experiment_metadata_yaml", "experiment_metadata.yaml")

    if not rrna_mapping_csv.exists():
        raise SystemExit(f"rRNA mapping CSV not found: {rrna_mapping_csv}")
    if not mag_annotation_dir.exists():
        raise SystemExit(f"Annotation directory not found: {mag_annotation_dir}")

    rrna_df = pd.read_csv(rrna_mapping_csv)

    experiment_meta: Dict = {}
    exp_meta_path = Path(experiment_metadata_yaml)
    if exp_meta_path.exists():
        with exp_meta_path.open() as fh:
            experiment_meta = yaml.safe_load(fh) or {}

    if meta_mapping_csv.exists():
        mapping_df = pd.read_csv(meta_mapping_csv)
        has_mapping = not mapping_df.empty
    else:
        mapping_df = pd.DataFrame()
        has_mapping = False

    # Aggregate abundances by MAG if mapping is available
    mag_abundance: Dict[str, Dict[str, float]] = {}
    abundance_cols: List[str] = []
    if has_mapping and "mag_id" in mapping_df.columns:
        abundance_cols = [
            c
            for c in mapping_df.columns
            if c
            not in {
                "metabarcoding_id",
                "metabarcoding_sequence",
                "rrna_uid",
                "mag_id",
                "sequence_header",
                "rrna_type",
                "identity",
            }
            and pd.api.types.is_numeric_dtype(mapping_df[c])
        ]
        if abundance_cols:
            grouped = mapping_df.groupby("mag_id")[abundance_cols].sum()
            mag_abundance = grouped.to_dict(orient="index")

    lines: List[str] = []
    # For community-level aggregation
    all_mag_ec: Dict[str, Dict[str, int]] = {}
    all_mag_cog: Dict[str, Dict[str, int]] = {}

    lines.append("# CoSMIC Composition Report")
    lines.append("")

    # MAG metadata section
    lines.append("## MAG Metadata")
    mag_ids = sorted(rrna_df["mag_id"].dropna().unique().tolist())
    if not mag_ids:
        lines.append("No MAGs found in rRNA mapping.")
    else:
        for mag_id in mag_ids:
            mag_rows = rrna_df[rrna_df["mag_id"] == mag_id]
            mag_path = mag_rows["mag_fasta_path"].iloc[0]
            n_rrna_16S = (mag_rows["rrna_type"] == "16S").sum()
            n_rrna_18S = (mag_rows["rrna_type"] == "18S").sum()

            ann_dir = mag_annotation_dir / mag_id
            summary_txt = ann_dir / f"{mag_id}.txt"
            summary_lines: List[str] = []
            if summary_txt.exists():
                with summary_txt.open() as fh:
                    summary_lines = [ln.strip() for ln in fh if ln.strip()]

            lines.append(f"### MAG {mag_id}")
            lines.append(f"- Assembly path: {mag_path}")
            lines.append(f"- rRNA features: {n_rrna_16S} × 16S, {n_rrna_18S} × 18S")
            if summary_lines:
                for ln in summary_lines:
                    lines.append(f"- {ln}")
            lines.append("")

    # Experiment metadata
    lines.append("## CoSMIC Experiment Metadata")
    if experiment_meta:
        for key, value in experiment_meta.items():
            lines.append(f"- {key}: {value}")
    else:
        lines.append("- (No experiment metadata file found; "
                     "add details to experiment_metadata.yaml)")
    lines.append("")

    # Species / MAG annotation + abundance
    lines.append("## Species / MAG Annotations and Relative Abundances")
    if not mag_ids:
        lines.append("No MAGs to summarize.")
    else:
        for mag_id in mag_ids:
            lines.append(f"### MAG {mag_id}")

            # Abundance info
            abund_info = mag_abundance.get(mag_id)
            if abund_info and abundance_cols:
                lines.append("- Relative abundance per sample (CoSMIC metabarcoding):")
                for col in abundance_cols:
                    val = abund_info.get(col, 0.0)
                    lines.append(f"  - {col}: {val:.4f}")
            else:
                lines.append("- Relative abundance: not available (no metabarcoding mapping file found).")

            # Functional annotation overview from Prokka TSV
            ann_dir = mag_annotation_dir / mag_id
            tsv_path = ann_dir / f"{mag_id}.tsv"
            if tsv_path.exists():
                tsv_df = pd.read_csv(tsv_path, sep="\t")
                cds_df = tsv_df[tsv_df["ftype"] == "CDS"]
                total_cds = len(cds_df)
                non_hyp_df = cds_df[~cds_df["product"].str.contains("hypothetical", case=False, na=False)]

                lines.append(f"- Annotated CDS count: {total_cds}")
                lines.append(f"- Non-hypothetical CDS count: {len(non_hyp_df)}")

                # Top products by frequency (brief overview)
                product_counts = (
                    non_hyp_df["product"]
                    .value_counts()
                    .head(10)
                    .to_dict()
                )
                if product_counts:
                    lines.append("- Most frequent annotated products:")
                    for prod, count in product_counts.items():
                        lines.append(f"  - {prod} (n={count})")
                else:
                    lines.append("- No non-hypothetical products with counts to summarize.")

                # EC numbers and COGs per MAG
                mag_ec_counts: Dict[str, int] = {}
                mag_cog_counts: Dict[str, int] = {}

                if "EC_number" in cds_df.columns:
                    for val in cds_df["EC_number"].dropna():
                        if not isinstance(val, str):
                            continue
                        for ec in [x.strip() for x in str(val).replace(";", ",").split(",")]:
                            if not ec or ec in {"-", "NA"}:
                                continue
                            mag_ec_counts[ec] = mag_ec_counts.get(ec, 0) + 1

                if "COG" in cds_df.columns:
                    for val in cds_df["COG"].dropna():
                        if not isinstance(val, str):
                            continue
                        for cog in [x.strip() for x in str(val).replace(";", ",").split(",")]:
                            if not cog or cog in {"-", "NA"}:
                                continue
                            mag_cog_counts[cog] = mag_cog_counts.get(cog, 0) + 1

                all_mag_ec[mag_id] = mag_ec_counts
                all_mag_cog[mag_id] = mag_cog_counts

                if mag_ec_counts:
                    lines.append("- Top EC numbers (by CDS count):")
                    for ec, count in sorted(
                        mag_ec_counts.items(), key=lambda kv: kv[1], reverse=True
                    )[:10]:
                        lines.append(f"  - {ec} (n={count})")
                else:
                    lines.append("- No EC annotations available for this MAG.")

                if mag_cog_counts:
                    lines.append("- Top COGs (by CDS count):")
                    for cog, count in sorted(
                        mag_cog_counts.items(), key=lambda kv: kv[1], reverse=True
                    )[:10]:
                        lines.append(f"  - {cog} (n={count})")
                else:
                    lines.append("- No COG annotations available for this MAG.")
            else:
                lines.append("- Annotation TSV not found; only presence/absence known.")

            lines.append("")

    # Community-level aggregation of ECs and COGs, weighted by MAG abundance
    lines.append("## Community-Level Functional Summary")
    if not mag_ids:
        lines.append("No MAGs to summarize at community level.")
    else:
        community_ec_counts: Dict[str, float] = {}
        community_cog_counts: Dict[str, float] = {}

        for mag_id in mag_ids:
            mag_ec = all_mag_ec.get(mag_id, {})
            mag_cog = all_mag_cog.get(mag_id, {})
            abund_info = mag_abundance.get(mag_id)
            if abund_info and abundance_cols:
                weight = sum(abund_info.get(col, 0.0) for col in abundance_cols)
            else:
                weight = 1.0

            for ec, count in mag_ec.items():
                community_ec_counts[ec] = community_ec_counts.get(ec, 0.0) + count * weight
            for cog, count in mag_cog.items():
                community_cog_counts[cog] = community_cog_counts.get(cog, 0.0) + count * weight

        if community_ec_counts:
            lines.append("- Top community EC numbers (abundance-weighted):")
            for ec, score in sorted(
                community_ec_counts.items(), key=lambda kv: kv[1], reverse=True
            )[:20]:
                lines.append(f"  - {ec} (weighted count ~ {score:.2f})")
        else:
            lines.append("- No EC annotations available at community level.")

        if community_cog_counts:
            lines.append("- Top community COGs (abundance-weighted):")
            for cog, score in sorted(
                community_cog_counts.items(), key=lambda kv: kv[1], reverse=True
            )[:20]:
                lines.append(f"  - {cog} (weighted count ~ {score:.2f})")
        else:
            lines.append("- No COG annotations available at community level.")

    lines.append("")

    # LLM query prompt
    lines.append("## LLM Query")
    lines.append(
        "Given the MAG metadata, CoSMIC experiment metadata, species/MAG annotations, "
        "and relative abundances described above, answer the following:\n"
    )
    lines.append(
        "1. Based on this community composition, what types of environments or "
        "substrates is this population likely engaging with?\n"
        "2. What biosynthetic and metabolic processes are likely active in this "
        "microbial community, and which MAGs/species appear to contribute "
        "to these functions?"
    )

    out_path = output_dir / (args.output or "cosmic_llm_report.md")
    with out_path.open("w") as fh:
        fh.write("\n".join(lines))

    print(f"Wrote LLM-ready report to {out_path}")


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
    run_parser.add_argument(
        "--output-dir",
        default=None,
        help="Base output directory for this run (default: current directory or config output_dir).",
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
    extract_parser.add_argument(
        "--output-dir",
        default=None,
        help="Base output directory for this extraction run (default: current directory or config output_dir).",
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
        "--output-dir",
        default=None,
        help=(
            "Base output directory for this annotation smoke-test "
            "(default: current directory or config output_dir)."
        ),
    )
    annotate_test_parser.add_argument(
        "--num-mags",
        type=int,
        default=5,
        help="How many MAGs to annotate (default: 5).",
    )

    report_parser = subparsers.add_parser(
        "report",
        help=(
            "Generate a markdown report summarizing MAGs, experiment metadata, "
            "and MAG annotations + abundances for LLM-based interpretation."
        ),
    )
    report_parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to YAML config file (default: config.yaml).",
    )
    report_parser.add_argument(
        "--rrna-mapping",
        default="barrnap_rrna_mapping.csv",
        help="Path to rRNA mapping CSV (default: barrnap_rrna_mapping.csv).",
    )
    report_parser.add_argument(
        "--annotation-dir",
        default=None,
        help="Directory containing Prokka annotation subfolders (overrides config).",
    )
    report_parser.add_argument(
        "--metabarcoding-mapping",
        default="metabarcoding_to_MAG_mapping.csv",
        help="Path to metabarcoding-to-MAG mapping CSV (if available).",
    )
    report_parser.add_argument(
        "--output",
        default="cosmic_llm_report.md",
        help="Report filename within the output directory (default: cosmic_llm_report.md).",
    )
    report_parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Base output directory for this report "
            "(default: current directory or config output_dir)."
        ),
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
    elif args.command == "report":
        cmd_report(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
