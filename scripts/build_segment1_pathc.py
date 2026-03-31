#!/usr/bin/env python3
"""
Path C: Build segment 1 tree from Path C genome set (~550).
Extract genome IDs, fetch FASTA, normalize headers, then augur align -> tree -> refine -> traits -> export.
Output: auspice/h3n2_segment1_pathc.json (segment 1 only; same Nextstrain format as main build).
Run tune_auspice_meta.py on that file after for color scales etc.
Does not modify Path A or Path B code or outputs.
"""
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
RESULTS_PATHC = PROJECT_ROOT / "results" / "segment_1_pathc"
AUSPICE_DIR = PROJECT_ROOT / "auspice"
CONFIG = PROJECT_ROOT / "config" / "config.yaml"

FASTA_PATHC = DATA_DIR / "segment_1_h3n2_pathc.clean.fasta"
META_PATHC = DATA_DIR / "metadata_segment1_pathc.tsv"
IDS_PATHC = DATA_DIR / "genome_ids_segment1_pathc.txt"
OUT_JSON = AUSPICE_DIR / "h3n2_segment1_pathc.json"


def run(cmd: list[str], check: bool = True) -> subprocess.CompletedProcess:
    print(" ", " ".join(cmd))
    return subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=check)


def main() -> int:
    print("=== Path C: Segment 1 build (~550 genomes) -> h3n2_segment1_pathc.json ===\n")

    if not IDS_PATHC.exists() or not META_PATHC.exists():
        print("1. Extracting genome IDs for Path C...")
        r = subprocess.run(
            [sys.executable, str(PROJECT_ROOT / "scripts" / "extract_genome_ids_from_path_c.py")],
            cwd=str(PROJECT_ROOT),
        )
        if r.returncode != 0:
            return 1
        if not IDS_PATHC.exists():
            print("Error: extract did not create genome_ids_segment1_pathc.txt")
            return 1
    else:
        print("1. genome_ids_segment1_pathc.txt and metadata exist; skipping extract.")

    if not (DATA_DIR / "segment_1_h3n2_pathc.fasta").exists():
        print("2. Fetching FASTA for Path C genome IDs...")
        run([sys.executable, "scripts/fetch_fasta_by_genome_id.py", "1", "--genome-list", str(IDS_PATHC), "--output-suffix", "_pathc", "--workers", "32"])
    else:
        print("2. segment_1_h3n2_pathc.fasta exists; skipping fetch.")

    print("3. Normalizing FASTA headers...")
    run([sys.executable, "scripts/ensure_fasta_headers_genome_id.py", "--suffix", "_pathc"])

    if not FASTA_PATHC.exists():
        print("Error: segment_1_h3n2_pathc.clean.fasta not found")
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

    RESULTS_PATHC.mkdir(parents=True, exist_ok=True)
    aligned = RESULTS_PATHC / "aligned.fasta"
    tree = RESULTS_PATHC / "tree.nwk"
    tree_refined = RESULTS_PATHC / "tree-refined.nwk"
    traits = RESULTS_PATHC / "traits_host.json"

    print("4. Aligning (MAFFT)...")
    run(["augur", "align", "--sequences", str(FASTA_PATHC), "--output", str(aligned), "--nthreads", str(nthreads)])

    print("5. Building tree (IQ-TREE)...")
    run(["augur", "tree", "--alignment", str(aligned), "--output", str(tree)])

    print("6. Refining tree (TreeTime)...")
    run([
        "augur", "refine", "--tree", str(tree), "--alignment", str(aligned), "--metadata", str(META_PATHC),
        "--output-tree", str(tree_refined), "--clock-rate", str(clock), "--coalescent", str(coalescent),
    ])

    print("7. Inferring traits (host)...")
    run(["augur", "traits", "--tree", str(tree_refined), "--metadata", str(META_PATHC), "--output", str(traits), "--columns", "host"])

    print("8. Exporting to Auspice JSON...")
    AUSPICE_DIR.mkdir(parents=True, exist_ok=True)
    run([
        "augur", "export", "v2", "--tree", str(tree_refined), "--metadata", str(META_PATHC),
        "--node-data", str(traits), "--output", str(OUT_JSON),
    ])

    print("9. Normalizing geo_resolutions...")
    run(["node", str(PROJECT_ROOT / "scripts" / "normalize_geo_resolutions.js")])

    print(f"\nDone. Output: {OUT_JSON}")
    print("Run: python3 scripts/tune_auspice_meta.py h3n2_segment1_pathc.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
