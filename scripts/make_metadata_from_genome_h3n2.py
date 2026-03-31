"""
Path A: Convert genome_metadata_H3N2.tsv (from p3-all-genomes) to Nextstrain format.
Sets strain = genome_id for exact matching with FASTA headers.
Deduplicates by genome_id (one row per strain per segment).
"""
import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from location_to_country import get_region_for_country

PROJECT_ROOT = Path(__file__).resolve().parent.parent
GENOME_META = PROJECT_ROOT / "genome_metadata_H3N2.tsv"
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
OUTPUT_ALL = DATA_DIR / "metadata_h3n2.tsv"


def normalize_date(s: str) -> str:
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


def main():
    if not GENOME_META.exists():
        print(f"Error: {GENOME_META} not found.")
        print('Run: p3-all-genomes --eq="genome_name,Influenza A virus" --eq=subtype,H3N2 \\')
        print("  --attr genome_id,genome_name,strain,collection_date,isolation_country,geographic_group,host_common_name,segment \\")
        print("  > genome_metadata_H3N2.tsv")
        return 1

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    with open(GENOME_META, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        raw = list(reader)
        # Normalize keys: genome.genome_id -> genome_id
        rows = []
        for row in raw:
            r = {}
            for k, v in row.items():
                key = k.split(".")[-1] if "." in k else k
                r[key] = v
            rows.append(r)

    def get(row, *keys):
        for k in keys:
            v = row.get(k, "").strip().strip('"')
            if v:
                return v
        return ""

    out_cols = ["strain", "strain_name", "virus", "date", "region", "country", "host", "segment", "clade"]
    by_segment = {}
    seen_all = {}  # genome_id -> first row (for metadata_h3n2.tsv)

    for row in rows:
        gid = get(row, "genome_id", "genome id")
        if not gid:
            continue
        strain_name = get(row, "strain", "genome_name", "genome name")
        date = normalize_date(get(row, "collection_date", "collection date"))
        country = get(row, "isolation_country", "isolation country")
        region = get(row, "geographic_group", "geographic group")
        host = get(row, "host_common_name", "host common name", "host_common_name")
        segment = get(row, "segment")
        clade = get(row, "clade")

        if country and not region:
            region = get_region_for_country(country) or ""

        out_row = {
            "strain": gid,
            "strain_name": strain_name,
            "virus": "flu_h3n2",
            "date": date,
            "region": region,
            "country": country,
            "host": host,
            "segment": segment,
            "clade": clade,
        }

        if gid not in seen_all:
            seen_all[gid] = out_row

        seg = segment.strip('"')
        if seg:
            by_segment.setdefault(seg, {})[gid] = out_row

    # Write metadata_h3n2.tsv (all, deduplicated)
    with open(OUTPUT_ALL, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(seen_all.values())
    print(f"Wrote {len(seen_all)} rows to {OUTPUT_ALL}")

    # Write per-segment
    for seg in sorted(by_segment.keys(), key=lambda x: int(x) if x.isdigit() else 0):
        path = DATA_DIR / f"metadata_segment{seg}_h3n2.tsv"
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=out_cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(by_segment[seg].values())
        print(f"Segment {seg}: {len(by_segment[seg])} rows -> {path}")

    return 0


if __name__ == "__main__":
    exit(main())
