#!/usr/bin/env python3
"""
Path C: Build segment N tree (N=1-8) from Path C genome set.
Extract genome IDs, fetch FASTA, normalize headers, then augur align -> tree -> refine -> traits -> export.
Output: auspice/h3n2_segmentN_pathc.json (one per segment).
Usage: python build_segment_pathc.py [--segment N]  (default: 1)
       python build_segment_pathc.py --all  (build segments 1-8)
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
AUSPICE_DIR = PROJECT_ROOT / "auspice"
CONFIG = PROJECT_ROOT / "config" / "config.yaml"


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    print(" ", " ".join(cmd))
    return subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=check)


def _count_fasta_headers(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def _resolve_pathc_input_fasta(seg: int) -> Path | None:
    """Path C input FASTA is Path B-style segment_N.clean.fasta (may live in two places)."""
    p1 = DATA_DIR / f"segment_{seg}.clean.fasta"
    if p1.exists():
        return p1
    p2 = PROJECT_ROOT / "H3N2_output" / f"segment_{seg}.clean.fasta"
    if p2.exists():
        return p2
    return None


def _read_lines(path: Path) -> list[str]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8", errors="replace") as f:
        return [ln.strip() for ln in f if ln.strip()]


def _clade_stats_from_metadata(meta_tsv: Path) -> dict[str, int]:
    """
    Compute unique clade/subclade counts from metadata TSV.
    We treat empty / 'unassigned' (case-insensitive) as missing.
    """
    stats = {
        "rows": 0,
        "unique_clades": 0,
        "unique_subclades": 0,
        "missing_clade": 0,
        "missing_subclade": 0,
    }
    if not meta_tsv.exists():
        return stats

    clades: set[str] = set()
    subclades: set[str] = set()
    missing_clade = 0
    missing_subclade = 0

    with open(meta_tsv, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            stats["rows"] += 1
            c = (row.get("clade") or "").strip().strip('"')
            sc = (row.get("subclade") or "").strip().strip('"')

            if not c or c.lower() == "unassigned":
                missing_clade += 1
            else:
                clades.add(c)

            if not sc or sc.lower() == "unassigned":
                missing_subclade += 1
            else:
                subclades.add(sc)

    stats["unique_clades"] = len(clades)
    stats["unique_subclades"] = len(subclades)
    stats["missing_clade"] = missing_clade
    stats["missing_subclade"] = missing_subclade
    return stats


def _print_summary_table(rows: list[dict]) -> None:
    if not rows:
        return

    headers = [
        "seg",
        "input_headers",
        "mapped_genome_ids",
        "metadata_rows",
        "clean_fasta_seqs",
        "unique_clades",
        "unique_subclades",
        "missing_clade",
        "missing_subclade",
    ]
    # Column widths
    widths = {h: len(h) for h in headers}
    for r in rows:
        for h in headers:
            widths[h] = max(widths[h], len(str(r.get(h, ""))))

    def fmt_row(r: dict) -> str:
        return "  ".join(str(r.get(h, "")).rjust(widths[h]) for h in headers)

    print("\n" + "=" * 60)
    print("Path C summary (metadata-based clade counts)")
    print("=" * 60)
    print("  " + "  ".join(h.rjust(widths[h]) for h in headers))
    print("  " + "  ".join("-" * widths[h] for h in headers))
    for r in rows:
        print("  " + fmt_row(r))
    print("=" * 60 + "\n")


def build_segment(seg: int) -> int:
    results_pathc = PROJECT_ROOT / "results" / f"segment_{seg}_pathc"
    fasta_pathc = DATA_DIR / f"segment_{seg}_h3n2_pathc.clean.fasta"
    meta_pathc = DATA_DIR / f"metadata_segment{seg}_pathc.tsv"
    ids_pathc = DATA_DIR / f"genome_ids_segment{seg}_pathc.txt"
    out_json = AUSPICE_DIR / f"h3n2_segment{seg}_pathc.json"

    if not ids_pathc.exists() or not meta_pathc.exists():
        print(f"1. Extracting genome IDs for Path C segment {seg}...")
        r = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "extract_genome_ids_from_path_c.py"), "--segment", str(seg)],
            cwd=str(PROJECT_ROOT),
        )
        if r.returncode != 0:
            return 1
        if not ids_pathc.exists():
            print(f"Error: extract did not create genome_ids_segment{seg}_pathc.txt")
            return 1
    else:
        print(f"1. genome_ids_segment{seg}_pathc.txt and metadata exist; skipping extract.")

    if not (DATA_DIR / f"segment_{seg}_h3n2_pathc.fasta").exists():
        print(f"2. Fetching FASTA for Path C segment {seg}...")
        run([sys.executable, "scripts/fetch_fasta_by_genome_id.py", str(seg), "--genome-list", str(ids_pathc), "--output-suffix", "_pathc", "--workers", "32"])
    else:
        print(f"2. segment_{seg}_h3n2_pathc.fasta exists; skipping fetch.")

    print("3. Normalizing FASTA headers...")
    run([sys.executable, "scripts/ensure_fasta_headers_genome_id.py", "--suffix", "_pathc", "--segment", str(seg)])

    if not fasta_pathc.exists():
        print(f"Error: segment_{seg}_h3n2_pathc.clean.fasta not found")
        return 1

    config = {}
    if CONFIG.exists():
        with open(CONFIG) as f:
            for line in f:
                m = re.match(r"(\w+):\s*(.+)", line.strip())
                if m:
                    config[m.group(1).strip()] = m.group(2).strip().strip('"')
    clock = config.get("clock_rate", "0.0008")
    coalescent = config.get("coalescent", "skyline")
    nthreads = config.get("align_nthreads", "auto")

    results_pathc.mkdir(parents=True, exist_ok=True)
    aligned = results_pathc / "aligned.fasta"
    tree = results_pathc / "tree.nwk"
    tree_refined = results_pathc / "tree-refined.nwk"
    traits = results_pathc / "traits_host.json"

    print("4. Aligning (MAFFT)...")
    run(["augur", "align", "--sequences", str(fasta_pathc), "--output", str(aligned), "--nthreads", str(nthreads)])

    print("5. Building tree (IQ-TREE)...")
    run(["augur", "tree", "--alignment", str(aligned), "--output", str(tree)])

    print("6. Refining tree (TreeTime)...")
    run([
        "augur", "refine", "--tree", str(tree), "--alignment", str(aligned), "--metadata", str(meta_pathc),
        "--output-tree", str(tree_refined), "--clock-rate", str(clock), "--coalescent", str(coalescent),
    ])

    print("7. Inferring traits (host)...")
    run(["augur", "traits", "--tree", str(tree_refined), "--metadata", str(meta_pathc), "--output", str(traits), "--columns", "host"])

    print("8. Exporting to Auspice JSON...")
    AUSPICE_DIR.mkdir(parents=True, exist_ok=True)
    run([
        "augur", "export", "v2", "--tree", str(tree_refined), "--metadata", str(meta_pathc),
        "--node-data", str(traits), "--output", str(out_json),
    ])

    print("9. Normalizing geo_resolutions...")
    run(["node", str(PROJECT_ROOT / "scripts" / "normalize_geo_resolutions.js")])

    print(f"\nDone. Output: {out_json}")
    print(f"Run: python scripts/tune_auspice_meta.py h3n2_segment{seg}_pathc.json")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description="Path C: Build segment tree(s)")
    ap.add_argument("--segment", type=int, default=1, help="Segment 1-8 (default 1)")
    ap.add_argument("--all", action="store_true", help="Build all segments 1-8")
    args = ap.parse_args()

    if args.all:
        segments = list(range(1, 9))
    else:
        if not (1 <= args.segment <= 8):
            print("Error: segment must be 1-8", file=sys.stderr)
            return 1
        segments = [args.segment]

    summary_rows: list[dict] = []

    for seg in segments:
        print(f"\n{'='*60}")
        print(f"=== Path C: Segment {seg} -> h3n2_segment{seg}_pathc.json ===")
        print("="*60)
        if build_segment(seg) != 0:
            return 1

        # Post-build verification + clade counts from metadata TSV
        inp = _resolve_pathc_input_fasta(seg)
        input_headers = _count_fasta_headers(inp) if inp else 0
        mapped_ids = len(_read_lines(DATA_DIR / f"genome_ids_segment{seg}_pathc.txt"))
        clean_seqs = _count_fasta_headers(DATA_DIR / f"segment_{seg}_h3n2_pathc.clean.fasta")
        clade_stats = _clade_stats_from_metadata(DATA_DIR / f"metadata_segment{seg}_pathc.tsv")

        summary_rows.append(
            {
                "seg": seg,
                "input_headers": input_headers,
                "mapped_genome_ids": mapped_ids,
                "metadata_rows": clade_stats.get("rows", 0),
                "clean_fasta_seqs": clean_seqs,
                "unique_clades": clade_stats.get("unique_clades", 0),
                "unique_subclades": clade_stats.get("unique_subclades", 0),
                "missing_clade": clade_stats.get("missing_clade", 0),
                "missing_subclade": clade_stats.get("missing_subclade", 0),
            }
        )

    if args.all:
        _print_summary_table(summary_rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
