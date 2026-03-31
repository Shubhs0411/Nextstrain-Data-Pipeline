#!/usr/bin/env python3
"""
Path A: Fetch FASTA for each segment by genome_id using p3-genome-fasta (parallel).

For 70k+ genomes per segment, the CLI is very slow (~8+ hours per segment). Prefer
BV-BRC website: upload genome_ids_segmentN.txt as a Genome List, then Export ->
Nucleotide Sequences. Save as segment_N_h3n2.fasta. Much faster.

This script supports --resume: skips IDs already present in the output file.

Usage:
  python3 scripts/fetch_fasta_by_genome_id.py 1 --limit 3   # test
  python3 scripts/fetch_fasta_by_genome_id.py 1 --resume   # resume segment 1
  python3 scripts/fetch_fasta_by_genome_id.py --all --workers 32
"""
import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"


def fetch_one(gid: str, timeout: int = 120) -> tuple[str, str | None]:
    """Run p3-genome-fasta for one genome_id. Returns (gid, stdout or None on failure)."""
    try:
        r = subprocess.run(
            ["p3-genome-fasta", "--contig", gid],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if r.returncode != 0:
            return (gid, None)
        return (gid, r.stdout)
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return (gid, None)


def parse_existing_fasta(path: Path) -> dict[str, str]:
    """Parse FASTA file; return dict genome_id -> full entry (header + sequence)."""
    if not path.exists():
        return {}
    out = {}
    current_id = None
    current_lines = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    out[current_id] = "".join(current_lines)
                header = line[1:].strip()
                first = header.split("|")[0].split()[0].strip()
                if ".con." in first:
                    current_id = first.split(".con.")[0]
                else:
                    current_id = first
                current_lines = [line]
            else:
                current_lines.append(line)
        if current_id is not None:
            out[current_id] = "".join(current_lines)
    return out


def fetch_fasta_for_segment(
    segment: int,
    limit: int | None = None,
    test_suffix: str = "",
    workers: int = 32,
    resume: bool = False,
    genome_list_path: Path | None = None,
    output_suffix: str = "",
) -> int:
    """Fetch FASTA for genome IDs in genome_ids_segment{N}.txt or custom list. Returns 0 on success."""
    suffix = output_suffix or test_suffix
    ids_file = genome_list_path or (DATA_DIR / f"genome_ids_segment{segment}.txt")
    out_file = DATA_DIR / f"segment_{segment}_h3n2{suffix}.fasta"

    if not ids_file.exists():
        print(f"Error: {ids_file} not found. Run write_genome_id_lists.py or extract_genome_ids_from_path_b.py first.")
        return 1

    ids = []
    with open(ids_file, encoding="utf-8") as f:
        for line in f:
            gid = line.strip()
            if gid and not gid.startswith("genome."):
                ids.append(gid)

    if limit is not None:
        ids = ids[:limit]

    existing = {}
    if resume and out_file.exists():
        existing = parse_existing_fasta(out_file)
        ids_to_fetch = [i for i in ids if i not in existing]
        print(f"Segment {segment}: resuming, {len(existing)} already in {out_file}, {len(ids_to_fetch)} to fetch")
    else:
        ids_to_fetch = ids
        if not resume and out_file.exists():
            print(f"Segment {segment}: {out_file} exists; use --resume to skip already-fetched IDs")

    if not ids_to_fetch:
        print(f"Segment {segment}: already complete ({len(ids)} genomes in {out_file})")
        return 0

    total = len(ids)
    print(f"Fetching FASTA for segment {segment} ({len(ids_to_fetch)} genomes, {workers} workers)...")

    results = dict(existing)
    done = 0
    report_every = max(1, min(500, len(ids_to_fetch) // 20))
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_gid = {executor.submit(fetch_one, gid): gid for gid in ids_to_fetch}
        for future in as_completed(future_to_gid):
            gid = future_to_gid[future]
            try:
                _, stdout = future.result()
                results[gid] = stdout
            except Exception:
                results[gid] = None
            done += 1
            if done % report_every == 0 or (limit and done <= 5):
                print(f"  {done}/{len(ids_to_fetch)} completed")

    with open(out_file, "w") as outf:
        for gid in ids:
            block = results.get(gid)
            if block:
                outf.write(block)
            else:
                print(f"  Warning: no data for {gid}", file=sys.stderr)

    ok = sum(1 for gid in ids if results.get(gid))
    print(f"Wrote {out_file} ({ok}/{len(ids)} genomes)")
    return 0


def main():
    ap = argparse.ArgumentParser(description="Fetch FASTA by genome_id for Path A (parallel, resumable)")
    ap.add_argument("segment", nargs="?", type=int, help="Segment number 1-8 (or use --all)")
    ap.add_argument("--all", action="store_true", help="Run for all segments 1-8")
    ap.add_argument("--limit", type=int, default=None, help="For testing: only fetch first N genome IDs")
    ap.add_argument("--test", action="store_true", help="Use _test.fasta suffix for trial run")
    ap.add_argument("--workers", type=int, default=32, help="Parallel workers (default 32)")
    ap.add_argument("--resume", action="store_true", help="Skip IDs already in output file (resume interrupted run)")
    ap.add_argument("--genome-list", type=Path, default=None, help="Use this file instead of genome_ids_segmentN.txt")
    ap.add_argument("--output-suffix", type=str, default="", help="Suffix for output file, e.g. _pathb -> segment_N_h3n2_pathb.fasta")
    args = ap.parse_args()

    if args.all:
        segments = list(range(1, 9))
    elif args.segment is not None:
        if not (1 <= args.segment <= 8):
            print("Error: segment must be 1-8", file=sys.stderr)
            return 1
        segments = [args.segment]
    else:
        ap.print_help()
        return 0

    suffix = "_test" if args.test else args.output_suffix
    genome_list = args.genome_list
    if genome_list and not genome_list.is_absolute():
        genome_list = PROJECT_ROOT / genome_list
    for seg in segments:
        if fetch_fasta_for_segment(
            seg,
            limit=args.limit,
            test_suffix="_test" if args.test else "",
            workers=args.workers,
            resume=args.resume,
            genome_list_path=genome_list,
            output_suffix=suffix,
        ) != 0:
            return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
