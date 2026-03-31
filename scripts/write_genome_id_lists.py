"""
Path A: Write genome_ids_segment1.txt ... genome_ids_segment8.txt
for use with BV-BRC to download FASTA by genome_id.
"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GENOME_META = PROJECT_ROOT / "genome_metadata_H3N2.tsv"
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"


def main():
    if not GENOME_META.exists():
        print(f"Error: {GENOME_META} not found. Run make_metadata_from_genome_h3n2.py first.")
        return 1

    by_segment = {}
    with open(GENOME_META, encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        # Support "genome.genome_id" or "genome_id"
        idx_gid = next((i for i, h in enumerate(header) if h.endswith("genome_id") or h == "genome_id"), 0)
        idx_seg = next((i for i, h in enumerate(header) if h.endswith("segment") or h == "segment"), -1)

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(idx_gid, idx_seg if idx_seg >= 0 else 0):
                continue
            gid = parts[idx_gid].strip().strip('"')
            seg = parts[idx_seg].strip().strip('"') if idx_seg >= 0 else "1"
            if gid and not gid.startswith("genome."):
                by_segment.setdefault(seg, set()).add(gid)

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for seg in sorted(by_segment.keys(), key=lambda x: int(x) if x.isdigit() else 0):
        path = DATA_DIR / f"genome_ids_segment{seg}.txt"
        with open(path, "w") as f:
            for gid in sorted(by_segment[seg]):
                f.write(gid + "\n")
        print(f"Segment {seg}: {len(by_segment[seg])} genome_ids -> {path}")

    return 0


if __name__ == "__main__":
    exit(main())
