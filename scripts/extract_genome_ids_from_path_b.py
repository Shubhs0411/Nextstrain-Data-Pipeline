#!/usr/bin/env python3
"""
Extract genome IDs for Path B strains (~550): from Path B FASTA headers, map strain name
to genome_id using BVBRC_genome.txt (primary) and genome_metadata_H3N2.tsv (fallback), then
write genome_ids_segment1_pathb.txt and metadata_segment1_pathb.tsv for the small build.
"""
import csv
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
# Path B clean FASTA (strain-name headers) - can be in data dir or output dir
PATHB_FASTA = DATA_DIR / "segment_1.clean.fasta"
GENOME_META = PROJECT_ROOT / "genome_metadata_H3N2.tsv"
BVBRC_GENOME = DATA_DIR / "BVBRC_genome.txt"
OUT_IDS = DATA_DIR / "genome_ids_segment1_pathb.txt"
OUT_META = DATA_DIR / "metadata_segment1_pathb.tsv"
PATH_A_META = DATA_DIR / "metadata_segment1_h3n2.tsv"


def load_pathb_strains(fasta_path: Path) -> set[str]:
    """Strain names from Path B FASTA headers (first token after '>')."""
    if not fasta_path.exists():
        return set()
    strains = set()
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split("|")[0].split()[0].strip()
                if name:
                    strains.add(name)
    return strains


def _normalize_strain(s: str) -> str:
    """Normalize for matching: strip quotes, collapse spaces/underscores."""
    s = (s or "").strip().strip('"')
    return s.replace(" ", "_")


def load_genome_meta_segment1() -> dict[str, str]:
    """From genome_metadata_H3N2.tsv segment=1: strain_name -> genome_id.
    Also adds normalized keys (spaces<->underscores, A/Human/ variant) so Path B headers match.
    """
    out = {}
    if not GENOME_META.exists():
        return out
    with open(GENOME_META, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = {k.split(".")[-1] if "." in k else k: v for k, v in row.items()}
            seg = (row.get("segment") or "").strip().strip('"')
            if seg != "1":
                continue
            gid = (row.get("genome_id") or "").strip().strip('"')
            strain = (row.get("strain") or "").strip().strip('"')
            if not gid or not strain or gid.startswith("genome."):
                continue
            out[strain] = gid
            # Normalized keys so Path B headers (often with underscores) match
            out[_normalize_strain(strain)] = gid
            out[strain.replace("_", " ")] = gid
            if strain.startswith("A/Human/") and "/" in strain[8:]:
                alt = "A/" + strain[8:]
                out[alt] = gid
                out[_normalize_strain(alt)] = gid
    return out


def _add_strain_variants(out: dict[str, str], strain: str, gid: str) -> None:
    """Add strain and all matching variants to lookup."""
    s = (strain or "").strip().strip('"')
    if not s or not gid:
        return
    out[s] = gid
    out[_normalize_strain(s)] = gid
    out[s.replace("_", " ")] = gid
    out[s.replace(" ", "_")] = gid
    # Case-insensitive: BVBRC uses lowercase host (swine, duck) vs Path B (Swine, Duck)
    out[s.lower()] = gid
    out[_normalize_strain(s.lower())] = gid
    # Collapse hyphen-space: "Helmern- IDT14829" -> "Helmern-IDT14829"
    collapsed = re.sub(r"-\s+", "-", s)
    if collapsed != s:
        out[collapsed] = gid
        out[_normalize_strain(collapsed)] = gid
        out[collapsed.lower()] = gid
    # mallard duck / ruddy shelduck -> Duck (Path B uses "Duck")
    for old, new in [("mallard duck", "Duck"), ("ruddy shelduck", "Duck")]:
        if old in s.lower():
            alt = re.sub(re.escape(old), new, s, flags=re.IGNORECASE)
            if alt != s:
                out[alt] = gid
                out[_normalize_strain(alt)] = gid
            break
    # Path B cleaning: A/Human/ -> A/
    if s.startswith("A/Human/") and "/" in s[8:]:
        alt = "A/" + s[8:]
        out[alt] = gid
        out[_normalize_strain(alt)] = gid
    # Strip (H3N2) suffix if present
    base = re.sub(r"\(H3N2\)\s*$", "", s, flags=re.IGNORECASE).strip()
    if base != s:
        out[base] = gid
        out[_normalize_strain(base)] = gid


def load_bvbrc_strain_to_genome_id(segment: str = "1") -> dict[str, str]:
    """From BVBRC_genome.txt: Strain -> Genome ID for given Segment.
    Adds normalized variants (spaces/underscores, A/Human/, etc.) for matching.
    """
    out = {}
    if not BVBRC_GENOME.exists():
        return out
    try:
        with open(BVBRC_GENOME, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                seg = (row.get("Segment") or row.get("segment") or "").strip().strip('"')
                if str(seg) != str(segment):
                    continue
                gid = (row.get("Genome ID") or row.get("genome_id") or "").strip().strip('"')
                strain = (row.get("Strain") or row.get("strain") or "").strip().strip('"')
                if gid and strain and not gid.startswith("genome."):
                    _add_strain_variants(out, strain, gid)
    except Exception as e:
        print(f"  Warning: BVBRC_genome.txt read error: {e}", file=sys.stderr)
    return out


def main() -> int:
    fasta_path = PATHB_FASTA
    if not fasta_path.exists():
        alt = PROJECT_ROOT / "H3N2_output" / "segment_1.clean.fasta"
        if alt.exists():
            fasta_path = alt
        else:
            print(f"Error: Path B FASTA not found at {PATHB_FASTA} or {alt}")
            return 1

    strains = load_pathb_strains(fasta_path)
    print(f"Path B segment 1 FASTA: {len(strains)} strain names")

    # BVBRC_genome.txt is the Path B metadata source - use it FIRST (best match for Path B strains)
    strain_to_gid = {}
    if BVBRC_GENOME.exists():
        strain_to_gid = load_bvbrc_strain_to_genome_id("1")
        print(f"BVBRC_genome.txt: {len(set(strain_to_gid.values()))} genome IDs for segment 1")
    # Merge in genome_metadata_H3N2.tsv for any additional matches
    meta_map = load_genome_meta_segment1()
    for k, v in meta_map.items():
        if k not in strain_to_gid:
            strain_to_gid[k] = v
    print(f"Combined lookup: {len(set(strain_to_gid.values()))} unique genome IDs")

    genome_ids = []
    for s in strains:
        gid = (
            strain_to_gid.get(s)
            or strain_to_gid.get(_normalize_strain(s))
            or strain_to_gid.get(s.replace(" ", "_"))
            or strain_to_gid.get(s.replace("_", " "))
            or strain_to_gid.get(s.lower())
            or strain_to_gid.get(_normalize_strain(s.lower()))
        )
        if not gid and s.startswith("A/Human/"):
            alt = "A/" + s[8:]
            gid = strain_to_gid.get(alt) or strain_to_gid.get(_normalize_strain(alt))
        if not gid:
            base = re.sub(r"\(H3N2\)\s*$", "", s, flags=re.IGNORECASE).strip()
            if base != s:
                gid = strain_to_gid.get(base) or strain_to_gid.get(_normalize_strain(base))
        if gid:
            genome_ids.append(gid)
        else:
            print(f"  No genome_id for strain: {s[:60]}...", file=sys.stderr)

    genome_ids = list(dict.fromkeys(genome_ids))  # dedup preserve order
    print(f"Mapped to {len(genome_ids)} genome IDs")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_IDS, "w") as f:
        for gid in genome_ids:
            f.write(gid + "\n")
    print(f"Wrote {OUT_IDS}")

    # Build metadata_segment1_pathb.tsv: Path A rows where available, else BVBRC rows
    out_cols = ["strain", "strain_name", "virus", "date", "region", "country", "host", "segment", "clade"]
    gid_set = set(genome_ids)
    gid_to_row: dict[str, dict] = {}

    # 1. From Path A metadata (strain=genome_id)
    if PATH_A_META.exists():
        with open(PATH_A_META, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                strain = (row.get("strain") or "").strip().strip('"')
                if strain in gid_set:
                    gid_to_row[strain] = {k: (row.get(k) or "").strip().strip('"') for k in out_cols}

    # 2. From BVBRC_genome.txt for genome_ids not in Path A
    missing = gid_set - set(gid_to_row)
    if missing and BVBRC_GENOME.exists():
        try:
            from location_to_country import get_region_for_country
        except ImportError:
            get_region_for_country = lambda c: ""

        def _norm_date(s: str) -> str:
            s = (s or "").strip().strip('"')
            if not s:
                return ""
            s = re.sub(r"T\d{2}:\d{2}:\d{2}.*$", "", s)
            if re.match(r"^\d{4}-\d{2}-\d{2}$", s):
                return s
            if re.match(r"^\d{4}-\d{2}$", s):
                return s + "-01"
            if re.match(r"^\d{4}$", s):
                return s + "-01-01"
            return s

        with open(BVBRC_GENOME, encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                gid = (row.get("Genome ID") or "").strip().strip('"')
                if gid not in missing:
                    continue
                seg = (row.get("Segment") or "").strip().strip('"')
                if str(seg) != "1":
                    continue
                strain_name = (row.get("Strain") or "").strip().strip('"')
                date = _norm_date(row.get("Collection Date", ""))
                country = (row.get("Isolation Country") or "").strip().strip('"')
                region = (row.get("Geographic Group") or "").strip().strip('"')
                if country and not region:
                    region = get_region_for_country(country) or ""
                host = (row.get("Host Common Name") or row.get("Host Name") or "").strip().strip('"')
                clade = (row.get("Clade") or "").strip().strip('"')
                gid_to_row[gid] = {
                    "strain": gid,
                    "strain_name": strain_name,
                    "virus": "flu_h3n2",
                    "date": date,
                    "region": region,
                    "country": country,
                    "host": host,
                    "segment": "1",
                    "clade": clade,
                }
                missing.discard(gid)

        if missing:
            print(f"  Warning: {len(missing)} genome IDs have no metadata (Path A or BVBRC)", file=sys.stderr)

    out_rows = [gid_to_row[gid] for gid in genome_ids if gid in gid_to_row]
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_META, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {OUT_META} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
