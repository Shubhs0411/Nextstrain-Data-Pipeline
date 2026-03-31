#!/usr/bin/env python3
"""
Path C: Extract genome IDs for the Path C strain set (~550 per segment).
Reads Path B FASTA (segment_N.clean.fasta) headers, maps strain name to genome_id
using BVBRC_genome.txt (primary) and genome_metadata_H3N2.tsv (fallback).
Writes genome_ids_segmentN_pathc.txt and metadata_segmentN_pathc.tsv.
Usage: python extract_genome_ids_from_path_c.py [--segment N]  (default: 1)
"""
import argparse
import csv
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
GENOME_META = PROJECT_ROOT / "genome_metadata_H3N2.tsv"
BVBRC_GENOME = DATA_DIR / "BVBRC_genome.txt"
NEW_CLADE_FILE = PROJECT_ROOT / "NEW_genome_w_clade.txt"


def _clean_fasta_strain(name: str) -> str:
    """Apply same cleaning as clean_fasta_headers so segment 2–8 headers match BVBRC."""
    s = (name or "").strip().strip('"')
    if not s:
        return s
    # Remove outer parentheses: (A/Human/...) -> A/Human/...
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    # Remove (H3N2) suffix
    s = re.sub(r"\(H3N2\)\s*$", "", s, flags=re.IGNORECASE).strip()
    # Normalize BV-BRC: A/Human/... -> A/...
    if s.startswith("A/Human/") and len(s) > 8:
        s = "A/" + s[8:]
    return s


def load_pathc_strains(fasta_path: Path) -> set[str]:
    """Strain names from FASTA headers (first token after '>'), cleaned for BVBRC matching."""
    if not fasta_path.exists():
        return set()
    strains = set()
    with open(fasta_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split("|")[0].split()[0].strip()
                if name:
                    cleaned = _clean_fasta_strain(name)
                    if cleaned:
                        strains.add(cleaned)
    return strains


def _normalize_strain(s: str) -> str:
    s = (s or "").strip().strip('"')
    return s.replace(" ", "_")


def _pubmed_url(pub: str) -> str:
    """Convert BV-BRC Publication field to a PubMed URL when possible."""
    s = (pub or "").strip().strip('"')
    if not s:
        return ""
    m = re.search(r"\d+", s)
    if m:
        return f"https://www.ncbi.nlm.nih.gov/pubmed/{m.group()}"
    return s


def load_new_genome_w_clade() -> dict[str, dict]:
    """Load NEW_genome_w_clade.txt: normalized strain -> {host,country,region,year,clade}."""
    out: dict[str, dict] = {}
    if not NEW_CLADE_FILE.exists():
        return out
    with open(NEW_CLADE_FILE, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            raw_id = (parts[0] or "").strip().strip('"')
            host = (parts[1] or "").strip().strip('"') if len(parts) > 1 else ""
            country_region = (parts[2] or "").strip().strip('"') if len(parts) > 2 else ""
            year = (parts[3] or "").strip().strip('"') if len(parts) > 3 else ""
            clade = (parts[4] or "").strip().strip('"') if len(parts) > 4 else ""

            if not raw_id:
                continue
            # Normalize strain_id: (A/Human/...) -> A/Human/... and spaces -> underscores, drop (H3N2)
            clean = raw_id
            if clean.startswith("(") and clean.endswith(")"):
                clean = clean[1:-1]
            clean = clean.replace(" ", "_").replace("(H3N2)", "").strip()

            # For now we don't try to resolve LOCATION_TO_COUNTRY; we just store raw country/region string.
            country_val = country_region.replace(" ", "_")
            region_val = ""

            out[clean] = {
                "host": host,
                "country": country_val,
                "region": region_val,
                "year": year,
                "clade": clade,
            }
            # Also store without host prefix for matching (A/Human/X -> A/X)
            if clean.startswith("A/Human/") and "/" in clean[8:]:
                alt = "A/" + clean[8:]
                if alt not in out:
                    out[alt] = out[clean].copy()
    return out


def load_genome_meta_segment(segment: str) -> dict[str, str]:
    out = {}
    if not GENOME_META.exists():
        return out
    with open(GENOME_META, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            row = {k.split(".")[-1] if "." in k else k: v for k, v in row.items()}
            seg = (row.get("segment") or "").strip().strip('"')
            if str(seg) != str(segment):
                continue
            gid = (row.get("genome_id") or "").strip().strip('"')
            strain = (row.get("strain") or "").strip().strip('"')
            if not gid or not strain or gid.startswith("genome."):
                continue
            out[strain] = gid
            out[_normalize_strain(strain)] = gid
            out[strain.replace("_", " ")] = gid
            if strain.startswith("A/Human/") and "/" in strain[8:]:
                alt = "A/" + strain[8:]
                out[alt] = gid
                out[_normalize_strain(alt)] = gid
    return out


def _add_strain_variants(out: dict[str, str], strain: str, gid: str) -> None:
    s = (strain or "").strip().strip('"')
    if not s or not gid:
        return
    out[s] = gid
    out[_normalize_strain(s)] = gid
    out[s.replace("_", " ")] = gid
    out[s.replace(" ", "_")] = gid
    out[s.lower()] = gid
    out[_normalize_strain(s.lower())] = gid
    collapsed = re.sub(r"-\s+", "-", s)
    if collapsed != s:
        out[collapsed] = gid
        out[_normalize_strain(collapsed)] = gid
        out[collapsed.lower()] = gid
    for old, new in [("mallard duck", "Duck"), ("ruddy shelduck", "Duck")]:
        if old in s.lower():
            alt = re.sub(re.escape(old), new, s, flags=re.IGNORECASE)
            if alt != s:
                out[alt] = gid
                out[_normalize_strain(alt)] = gid
            break
    if s.startswith("A/Human/") and "/" in s[8:]:
        alt = "A/" + s[8:]
        out[alt] = gid
        out[_normalize_strain(alt)] = gid
    base = re.sub(r"\(H3N2\)\s*$", "", s, flags=re.IGNORECASE).strip()
    if base != s:
        out[base] = gid
        out[_normalize_strain(base)] = gid


def load_bvbrc_strain_to_genome_id(segment: str = "1") -> dict[str, str]:
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
    ap = argparse.ArgumentParser(description="Path C: Extract genome IDs from segment FASTA")
    ap.add_argument("--segment", type=int, default=1, help="Segment 1-8 (default 1)")
    args = ap.parse_args()
    seg = args.segment
    if not (1 <= seg <= 8):
        print("Error: segment must be 1-8", file=sys.stderr)
        return 1

    fasta_path = DATA_DIR / f"segment_{seg}.clean.fasta"
    if not fasta_path.exists():
        alt = PROJECT_ROOT / "H3N2_output" / f"segment_{seg}.clean.fasta"
        if alt.exists():
            fasta_path = alt
        else:
            print(f"Error: Path C input FASTA not found at {fasta_path} or {alt}")
            return 1

    out_ids = DATA_DIR / f"genome_ids_segment{seg}_pathc.txt"
    out_meta = DATA_DIR / f"metadata_segment{seg}_pathc.tsv"
    path_a_meta = DATA_DIR / f"metadata_segment{seg}_h3n2.tsv"

    strains = load_pathc_strains(fasta_path)
    print(f"Path C segment {seg} FASTA: {len(strains)} strain names")

    # Load clade overrides from NEW_genome_w_clade.txt (strain-based)
    new_clade = load_new_genome_w_clade()

    strain_to_gid = {}
    if BVBRC_GENOME.exists():
        strain_to_gid = load_bvbrc_strain_to_genome_id(str(seg))
        print(f"BVBRC_genome.txt: {len(set(strain_to_gid.values()))} genome IDs for segment {seg}")
    meta_map = load_genome_meta_segment(str(seg))
    for k, v in meta_map.items():
        if k not in strain_to_gid:
            strain_to_gid[k] = v
    print(f"Combined lookup: {len(set(strain_to_gid.values()))} unique genome IDs")
    if not strain_to_gid:
        print("Error: no strain->genome_id mapping found", file=sys.stderr)
        return 1

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

    genome_ids = list(dict.fromkeys(genome_ids))
    print(f"Mapped to {len(genome_ids)} unique genome IDs (from {len(strains)} input strains)")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(out_ids, "w") as f:
        for gid in genome_ids:
            f.write(gid + "\n")
    print(f"Wrote {out_ids}")

    # Keep standard Nextstrain columns plus extra BV-BRC fields we can show in Auspice
    out_cols = [
        "strain", "strain_name", "virus", "date",
        "region", "country", "division", "city",
        "host", "segment", "clade", "subclade", "lineage",
        "accession", "authors", "paper_url", "journal", "isolation_source",
    ]
    gid_set = set(genome_ids)
    gid_to_row: dict[str, dict] = {}

    if path_a_meta.exists():
        with open(path_a_meta, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                strain = (row.get("strain") or "").strip().strip('"')
                if strain in gid_set:
                    gid_to_row[strain] = {k: (row.get(k) or "").strip().strip('"') for k in out_cols}

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
                row_seg = (row.get("Segment") or "").strip().strip('"')
                if str(row_seg) != str(seg):
                    continue
                strain_name = (row.get("Strain") or "").strip().strip('"')
                date = _norm_date(row.get("Collection Date", ""))
                country = (row.get("Isolation Country") or "").strip().strip('"')
                region = (row.get("Geographic Group") or "").strip().strip('"')
                if country and not region:
                    region = get_region_for_country(country) or ""
                host = (row.get("Host Common Name") or row.get("Host Name") or "").strip().strip('"')
                clade = (row.get("Clade") or "").strip().strip('"')
                subclade = (row.get("Subclade") or "").strip().strip('"')
                lineage = (row.get("Lineage") or "").strip().strip('"')
                division = (row.get("State/Province") or "").strip().strip('"')
                city = (row.get("City") or "").strip().strip('"')
                isolation_source = (row.get("Isolation Source") or "").strip().strip('"')
                authors = (row.get("Authors") or "").strip().strip('"')
                pub = (row.get("Publication") or "").strip().strip('"')
                accession = (row.get("GenBank Accessions") or "").strip().strip('"')
                paper_url = _pubmed_url(pub)
                journal = ""  # BV-BRC Publication field is usually PMID, so journal is often not provided

                # Overlay clade / basic fields from NEW_genome_w_clade.txt when available
                if strain_name:
                    strain_clean = strain_name.replace(" ", "_").replace("(H3N2)", "")
                    nc = new_clade.get(strain_clean)
                    if nc:
                        if nc.get("clade") and nc["clade"].lower() != "unassigned":
                            clade = nc["clade"]
                        if nc.get("country") and not country:
                            country = nc["country"]
                        if nc.get("region") and not region:
                            region = nc["region"]
                        if nc.get("host") and not host:
                            host = nc["host"]
                        if nc.get("year") and not date:
                            y = nc["year"]
                            date = f"{y}-01-01" if len(y) == 4 else y
                gid_to_row[gid] = {
                    "strain": gid,
                    "strain_name": strain_name,
                    "virus": "flu_h3n2",
                    "date": date,
                    "region": region,
                    "country": country,
                    "division": division,
                    "city": city,
                    "host": host,
                    "segment": str(seg),
                    "clade": clade,
                    "subclade": subclade,
                    "lineage": lineage,
                    "accession": accession,
                    "authors": authors,
                    "paper_url": paper_url,
                    "journal": journal,
                    "isolation_source": isolation_source,
                }
                missing.discard(gid)

        if missing:
            print(f"  Warning: {len(missing)} genome IDs have no metadata", file=sys.stderr)

    # Build metadata rows in the same order as genome_ids, then sanity‑check coverage.
    out_rows = [gid_to_row[gid] for gid in genome_ids if gid in gid_to_row]
    if len(out_rows) != len(genome_ids):
        missing_meta = [gid for gid in genome_ids if gid not in gid_to_row]
        print(
            f"  Warning: {len(missing_meta)} genome IDs have no metadata rows "
            f"(metadata rows={len(out_rows)} vs genome_ids={len(genome_ids)})",
            file=sys.stderr,
        )
        if missing_meta:
            preview = ", ".join(missing_meta[:5])
            if len(missing_meta) > 5:
                preview += ", ..."
            print(f"    Example missing genome_ids: {preview}", file=sys.stderr)
    with open(out_meta, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {out_meta} ({len(out_rows)} rows; genome_ids={len(genome_ids)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
