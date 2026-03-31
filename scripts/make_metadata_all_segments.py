"""
Filter metadata by segment and deduplicate by strain.
BV-BRC has 8 rows per strain (one per segment); output one row per strain per segment file.
"""
import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "H3N2_output"
METADATA_FILE = OUTPUT_DIR / "metadata.tsv"
SEGMENTS = list(range(1, 9))


def main():
    if not METADATA_FILE.exists():
        print(f"Error: {METADATA_FILE} not found. Run make_metadata.py first.")
        return 1

    with open(METADATA_FILE, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames = reader.fieldnames or list(rows[0].keys()) if rows else []

    print(f"Loaded {len(rows)} rows from {METADATA_FILE}")

    for seg in SEGMENTS:
        seg_str = str(seg)
        seg_rows = [r for r in rows if r.get("segment", "").strip('"') == seg_str]
        # Deduplicate by strain (keep first occurrence)
        seen = set()
        unique = []
        for r in seg_rows:
            strain = r.get("strain", "").strip()
            if strain and strain not in seen:
                seen.add(strain)
                unique.append(r)

        out_path = OUTPUT_DIR / f"metadata_segment{seg}.tsv"
        with open(out_path, "w", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(unique)

        print(f"Segment {seg}: {len(unique)} unique strains -> {out_path}")

    return 0


if __name__ == "__main__":
    exit(main())
