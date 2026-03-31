"""
Path A: Normalize FASTA headers to genome_id and deduplicate.
Reads segment_N_h3n2.fasta (or segment_N.fasta), rewrites headers to genome_id.
Keeps only the first occurrence of each genome_id (dedup) so augur align succeeds.
Use --suffix _pathb to process segment_1_h3n2_pathb.fasta -> segment_1_h3n2_pathb.clean.fasta.
"""
import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
SEGMENTS = list(range(1, 9))


def load_valid_genome_ids(segment: int, suffix: str = "") -> set:
    if suffix:
        # Path B small build uses metadata_segment1_pathb.tsv (no _h3n2 in name)
        path = DATA_DIR / f"metadata_segment{segment}{suffix}.tsv"
    else:
        path = DATA_DIR / f"metadata_segment{segment}_h3n2.tsv"
    if not path.exists():
        path = DATA_DIR / f"metadata_segment{segment}_h3n2{suffix}.tsv" if suffix else DATA_DIR / f"metadata_segment{segment}.tsv"
    if not path.exists():
        return set()
    ids = set()
    with open(path, encoding="utf-8") as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                ids.add(parts[0].strip().strip('"'))
    return ids


def header_to_genome_id(h: str, valid: set) -> str:
    """Extract genome_id from header."""
    tok = h.split("|")[0].split()[0].strip()
    if tok in valid:
        return tok
    if ".con." in tok:
        gid = tok.split(".con.")[0]
        if gid in valid:
            return gid
    return tok


def main():
    ap = argparse.ArgumentParser(description="Normalize FASTA headers to genome_id")
    ap.add_argument("--suffix", type=str, default="", help="e.g. _pathb for segment_N_h3n2_pathb.fasta")
    ap.add_argument("--segment", type=int, default=None, help="With --suffix: process only this segment (1-8)")
    args = ap.parse_args()
    suffix = (args.suffix or "").strip()
    if suffix and args.segment is not None:
        segments = [args.segment] if 1 <= args.segment <= 8 else [1]
    else:
        segments = [1] if suffix else SEGMENTS
    for seg in segments:
        valid = load_valid_genome_ids(seg, suffix)
        if not valid:
            print(f"  Skip segment {seg}: no metadata (expected metadata_segment{seg}{suffix or '_h3n2'}.tsv)")
            continue

        names = [f"segment_{seg}_h3n2{suffix}.fasta", f"segment_{seg}_h3n2.fasta", f"segment_{seg}.fasta"]
        for name in names:
            inp = DATA_DIR / name
            if not inp.exists():
                continue
            out = DATA_DIR / name.replace(".fasta", ".clean.fasta")
            seen = set()
            count = 0
            skipped = 0
            current_id = None
            current_lines = []

            with open(inp, encoding="utf-8") as fin, open(out, "w", encoding="utf-8") as fout:
                for line in fin:
                    if line.startswith(">"):
                        # Write previous record if first occurrence of that genome_id
                        if current_id is not None:
                            if current_id not in seen:
                                seen.add(current_id)
                                for l in current_lines:
                                    fout.write(l)
                                count += 1
                            else:
                                skipped += 1

                        h = line[1:].strip()
                        gid = header_to_genome_id(h, valid)
                        current_id = gid
                        current_lines = [">" + gid + "\n"]
                    else:
                        current_lines.append(line)

                # Last record
                if current_id is not None:
                    if current_id not in seen:
                        seen.add(current_id)
                        for l in current_lines:
                            fout.write(l)
                        count += 1
                    else:
                        skipped += 1

            if count > 0 or skipped > 0:
                dup_msg = f" ({skipped} duplicates skipped)" if skipped else ""
                print(f"Segment {seg}: {count} sequences -> {out}{dup_msg}")
                # Sanity check: how many genome_ids had metadata vs sequences written
                if len(valid) != count:
                    print(
                        f"  Warning: metadata genome_ids={len(valid)} but clean FASTA has {count} sequences.",
                        file=sys.stderr,
                    )
            break
    return 0


if __name__ == "__main__":
    exit(main())
