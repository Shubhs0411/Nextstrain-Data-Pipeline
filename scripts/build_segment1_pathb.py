#!/usr/bin/env python3
"""
Build segment 1 tree from Path B genome set (~500): extract genome IDs, fetch FASTA,
normalize headers, then run augur align -> tree -> refine -> traits -> export.
Output: auspice/h3n2_segment1_pathb.json (same Nextstrain format as main build).
Run tune_auspice_meta.py on that file after for color scales etc.
"""
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
RESULTS_PATHB = PROJECT_ROOT / "results" / "segment_1_pathb"
AUSPICE_DIR = PROJECT_ROOT / "auspice"
CONFIG = PROJECT_ROOT / "config" / "config.yaml"

FASTA_PATHB = DATA_DIR / "segment_1_h3n2_pathb.clean.fasta"
META_PATHB = DATA_DIR / "metadata_segment1_pathb.tsv"
IDS_PATHB = DATA_DIR / "genome_ids_segment1_pathb.txt"
OUT_JSON = AUSPICE_DIR / "h3n2_segment1_pathb.json"


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    print(" ", " ".join(cmd))
    return subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=check)


def main() -> int:
    print("=== Segment 1 Path B build (~500 genomes) -> h3n2_segment1_pathb.json ===\n")

    # 1. Extract genome IDs from Path B (if not already done)
    if not IDS_PATHB.exists() or not META_PATHB.exists():
        print("1. Extracting genome IDs from Path B...")
        r = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "extract_genome_ids_from_path_b.py")],
            cwd=str(PROJECT_ROOT),
        )
        if r.returncode != 0:
            return 1
        if not IDS_PATHB.exists():
            print("Error: extract did not create genome_ids_segment1_pathb.txt")
            return 1
    else:
        print("1. genome_ids_segment1_pathb.txt and metadata subset exist; skipping extract.")

    # 2. Fetch FASTA for those genome IDs
    if not (DATA_DIR / "segment_1_h3n2_pathb.fasta").exists():
        print("2. Fetching FASTA for Path B genome IDs...")
        run([sys.executable, "scripts/fetch_fasta_by_genome_id.py", "1", "--genome-list", str(IDS_PATHB), "--output-suffix", "_pathb", "--workers", "32"])
    else:
        print("2. segment_1_h3n2_pathb.fasta exists; skipping fetch.")

    # 3. Normalize FASTA headers to genome_id
    print("3. Normalizing FASTA headers...")
    run([sys.executable, "scripts/ensure_fasta_headers_genome_id.py", "--suffix", "_pathb"])

    if not FASTA_PATHB.exists():
        print("Error: segment_1_h3n2_pathb.clean.fasta not found")
        return 1

    # Load config (simple key: value)
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

    RESULTS_PATHB.mkdir(parents=True, exist_ok=True)
    aligned = RESULTS_PATHB / "aligned.fasta"
    tree = RESULTS_PATHB / "tree.nwk"
    tree_refined = RESULTS_PATHB / "tree-refined.nwk"
    traits = RESULTS_PATHB / "traits_host.json"

    # 4. Align
    print("4. Aligning (MAFFT)...")
    run(["augur", "align", "--sequences", str(FASTA_PATHB), "--output", str(aligned), "--nthreads", str(nthreads)])

    # 5. Tree
    print("5. Building tree (IQ-TREE)...")
    run(["augur", "tree", "--alignment", str(aligned), "--output", str(tree)])

    # 6. Refine
    print("6. Refining tree (TreeTime)...")
    run([
        "augur", "refine", "--tree", str(tree), "--alignment", str(aligned), "--metadata", str(META_PATHB),
        "--output-tree", str(tree_refined), "--clock-rate", str(clock), "--coalescent", str(coalescent),
    ])

    # 7. Traits
    print("7. Inferring traits (host)...")
    run(["augur", "traits", "--tree", str(tree_refined), "--metadata", str(META_PATHB), "--output", str(traits), "--columns", "host"])

    # 8. Export
    print("8. Exporting to Auspice JSON...")
    AUSPICE_DIR.mkdir(parents=True, exist_ok=True)
    run([
        "augur", "export", "v2", "--tree", str(tree_refined), "--metadata", str(META_PATHB),
        "--node-data", str(traits), "--output", str(OUT_JSON),
    ])

    print("9. Normalizing geo_resolutions...")
    run(["node", str(PROJECT_ROOT / "scripts" / "normalize_geo_resolutions.js")])

    print(f"\nDone. Output: {OUT_JSON}")
    print("Run: python3 scripts/tune_auspice_meta.py h3n2_segment1_pathb.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
