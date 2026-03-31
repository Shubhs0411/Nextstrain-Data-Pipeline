#!/usr/bin/env python3
"""
Verify Path C for segments 1-8: compare FASTA sequence counts vs mapped genome IDs and outputs.
Run after build_segment_pathc.py --all.
"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
OUTPUT_DIR = PROJECT_ROOT / "H3N2_output"
AUSPICE_DIR = PROJECT_ROOT / "auspice"


def count_fasta_headers(path: Path) -> int:
    if not path.exists():
        return -1
    n = 0
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def count_lines(path: Path) -> int:
    if not path.exists():
        return -1
    with open(path, encoding="utf-8") as f:
        return sum(1 for _ in f)


def main() -> int:
    print("Path C verification (segments 1-8)\n")
    print(f"{'Seg':<4} {'FASTA_in':<10} {'genome_ids':<12} {'metadata':<10} {'fetched':<10} {'clean':<10} {'JSON':<8}")
    print("-" * 70)

    for seg in range(1, 9):
        # Input FASTA (Path C reads from .clean.fasta or fallback)
        fasta_in = DATA_DIR / f"segment_{seg}.clean.fasta"
        if not fasta_in.exists():
            fasta_in = OUTPUT_DIR / f"segment_{seg}.clean.fasta"
        n_fasta = count_fasta_headers(fasta_in)

        ids_file = DATA_DIR / f"genome_ids_segment{seg}_pathc.txt"
        meta_file = DATA_DIR / f"metadata_segment{seg}_pathc.tsv"
        fetched = DATA_DIR / f"segment_{seg}_h3n2_pathc.fasta"
        clean_fasta = DATA_DIR / f"segment_{seg}_h3n2_pathc.clean.fasta"
        json_file = AUSPICE_DIR / f"h3n2_segment{seg}_pathc.json"

        n_ids = count_lines(ids_file)
        n_meta = (count_lines(meta_file) - 1) if meta_file.exists() else -1  # minus header
        if n_meta < 0 and meta_file.exists():
            n_meta = 0
        n_fetched = count_fasta_headers(fetched)
        n_clean = count_fasta_headers(clean_fasta)
        json_ok = "yes" if json_file.exists() else "no"

        print(f"{seg:<4} {n_fasta:<10} {n_ids:<12} {n_meta:<10} {n_fetched:<10} {n_clean:<10} {json_ok:<8}")

    print("-" * 70)
    print("\nColumns:")
    print("  FASTA_in   = sequences in segment_N.clean.fasta (source)")
    print("  genome_ids = lines in genome_ids_segmentN_pathc.txt (mapped)")
    print("  metadata   = data rows in metadata_segmentN_pathc.tsv")
    print("  fetched    = sequences in segment_N_h3n2_pathc.fasta")
    print("  clean      = sequences in segment_N_h3n2_pathc.clean.fasta")
    print("  JSON       = h3n2_segmentN_pathc.json exists")
    return 0


if __name__ == "__main__":
    sys.exit(main())
