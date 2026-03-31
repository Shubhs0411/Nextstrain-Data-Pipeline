#!/usr/bin/env python3
"""
Path C optional: overlay Nextclade clade / subclade onto Path C metadata.

Workflow:
  1) Run Path C build so that segment_4_h3n2_pathc.fasta exists.
  2) In your nextclade conda env, run (example):

       nextclade dataset get \
         --name flu_h3n2_ha \
         --output-dir nextclade_dataset_h3n2_ha

       nextclade run \
         --input-dataset nextclade_dataset_h3n2_ha \
         --output-all=H3N2_output/ \
         --output-basename nextclade_pathc_segment4 \
         H3N2_DATA/H3N2_DATA/segment_4_h3n2_pathc.fasta

     This produces H3N2_output/nextclade_pathc_segment4.tsv.

  3) In your bvbrc env (H3N2 project root), run:

       python scripts/add_nextclade_to_pathc.py \
         --tsv H3N2_output/nextclade_pathc_segment4.tsv

     This script:
       - Parses the Nextclade TSV.
       - Extracts genome_id from seqName (e.g. 518978.3 from 518978.3.con.0001 segment).
       - Updates clade / subclade columns in:
           H3N2_DATA/H3N2_DATA/metadata_segmentN_pathc.tsv  (N=1-8).

  4) Re-run tune_auspice_meta so clades appear in JSON:

       python scripts/tune_all_auspice_meta.py
     or
       python scripts/tune_auspice_meta.py h3n2_segmentN_pathc.json
"""

import argparse
import csv
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
OUTPUT_DIR = PROJECT_ROOT / "H3N2_output"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Overlay Nextclade clade/subclade onto Path C metadata_segmentN_pathc.tsv"
    )
    ap.add_argument(
        "--tsv",
        type=str,
        default=str(OUTPUT_DIR / "nextclade_pathc_segment4.tsv"),
        help="Path to Nextclade TSV output (default: H3N2_output/nextclade_pathc_segment4.tsv)",
    )
    ap.add_argument(
        "--segments",
        type=str,
        default="1-8",
        help="Segments to update, e.g. '1-8' or '1,2,4' (default: 1-8)",
    )
    return ap.parse_args()


def parse_segment_list(spec: str) -> list[int]:
    spec = (spec or "").strip()
    if not spec:
        return list(range(1, 9))
    if "-" in spec:
        start, end = spec.split("-", 1)
        try:
            s = int(start)
            e = int(end)
        except ValueError:
            return list(range(1, 9))
        s = max(1, s)
        e = min(8, e)
        if s > e:
            s, e = e, s
        return list(range(s, e + 1))
    parts = []
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            n = int(tok)
        except ValueError:
            continue
        if 1 <= n <= 8:
            parts.append(n)
    return sorted(set(parts)) or list(range(1, 9))


def header_to_genome_id(seq_name: str) -> str:
    """
    Extract genome_id from Nextclade seqName.
    Examples:
      '518978.3.con.0001 segment' -> '518978.3'
      '1184565.3.con.0001'        -> '1184565.3'
    """
    s = (seq_name or "").strip().strip('"')
    if not s:
        return ""
    tok = s.split()[0]
    if ".con." in tok:
        return tok.split(".con.")[0]
    return tok


def load_nextclade_map(tsv_path: Path) -> dict[str, dict[str, str]]:
    """
    Load Nextclade TSV. Returns:
      genome_id -> {"clade": ..., "subclade": ...}

    We map Nextclade columns to Path C fields as:
      - Path C clade <- legacy-clade (preferred), then short-clade, then clade.
      - Path C subclade <- clade (preferred), then subclade, then proposedSubclade.

    Rationale:
      For H3N2 HA, Nextclade often uses clade labels like J.2 while
      legacy-clade/short-clade carry 3C.2a-style labels. This keeps:
        clade ~ 3C lineage
        subclade ~ J.x lineage
    """
    def _clean(v: str) -> str:
        s = (v or "").strip().strip('"')
        if not s or s.lower() == "unassigned":
            return ""
        return s

    m: dict[str, dict[str, str]] = {}
    if not tsv_path.exists():
        print(f"Error: Nextclade TSV not found at {tsv_path}", file=sys.stderr)
        return m

    with open(tsv_path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq_name = row.get("seqName") or row.get("seqname") or ""
            gid = header_to_genome_id(seq_name)
            if not gid:
                continue
            nc_clade = _clean(row.get("clade") or "")
            nc_sub = _clean(row.get("subclade") or "")
            nc_prop = _clean(row.get("proposedSubclade") or "")
            nc_short = _clean(row.get("short-clade") or "")
            nc_legacy = _clean(row.get("legacy-clade") or "")

            clade = nc_legacy or nc_short or nc_clade
            sub = nc_clade or nc_sub or nc_prop

            if not clade and not sub:
                # Nothing new to contribute
                continue

            existing = m.get(gid, {})
            # Prefer non-empty clade/subclade, but don't overwrite with empty
            if clade:
                existing["clade"] = clade
            if sub:
                existing["subclade"] = sub
            m[gid] = existing

    print(f"Loaded Nextclade annotations for {len(m)} genome IDs from {tsv_path}")
    return m


def update_metadata_for_segment(seg: int, nc_map: dict[str, dict[str, str]]) -> None:
    meta_path = DATA_DIR / f"metadata_segment{seg}_pathc.tsv"
    if not meta_path.exists():
        print(f"  Segment {seg}: {meta_path.name} not found, skipping.")
        return

    rows = []
    updated = 0
    total = 0
    unique_clades: set[str] = set()
    unique_subclades: set[str] = set()

    with open(meta_path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        # Ensure clade / subclade columns exist
        if "clade" not in fieldnames:
            fieldnames.append("clade")
        if "subclade" not in fieldnames:
            fieldnames.append("subclade")

        for row in reader:
            total += 1
            gid = (row.get("strain") or "").strip().strip('"')
            info = nc_map.get(gid)
            if info:
                # Overlay clade / subclade from Nextclade if present
                clade = info.get("clade", "").strip()
                sub = info.get("subclade", "").strip()
                if clade:
                    row["clade"] = clade
                if sub:
                    row["subclade"] = sub
                updated += 1

            c = (row.get("clade") or "").strip().strip('"')
            sc = (row.get("subclade") or "").strip().strip('"')
            if c and c.lower() != "unassigned":
                unique_clades.add(c)
            if sc and sc.lower() != "unassigned":
                unique_subclades.add(sc)

            rows.append(row)

    # Write back in-place
    with open(meta_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"  Segment {seg}: updated {updated}/{total} rows "
        f"with Nextclade clade/subclade "
        f"(unique clades={len(unique_clades)}, unique subclades={len(unique_subclades)})"
    )


def main() -> int:
    args = parse_args()
    tsv_path = Path(args.tsv)
    segments = parse_segment_list(args.segments)

    nc_map = load_nextclade_map(tsv_path)
    if not nc_map:
        return 1

    print(f"Updating Path C metadata_segmentN_pathc.tsv for segments: {segments}")
    for seg in segments:
        update_metadata_for_segment(seg, nc_map)

    print("\nDone. Now re-run tune_auspice_meta so clades appear in JSON, e.g.:")
    print("  python scripts/tune_all_auspice_meta.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())

