"""
Path B optional: Fetch metadata for FASTA-only strains from BV-BRC.
Reads genome_metadata_H3N2.tsv (from p3-all-genomes) or runs p3-all-genomes,
matches by strain name, writes bvbrc_enriched.tsv for tune_auspice_meta.
"""
import csv
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from location_to_country import get_region_for_country, LOCATION_TO_COUNTRY

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "H3N2_output"
GENOME_META = PROJECT_ROOT / "genome_metadata_H3N2.tsv"
REPORT_FILE = OUTPUT_DIR / "fasta_metadata_comparison_report.txt"
ENRICHED_FILE = OUTPUT_DIR / "bvbrc_enriched.tsv"

KNOWN_HOSTS = {
    "Human", "Swine", "Duck", "Chicken", "Canine", "Turkey", "Environment",
    "Feline", "Mallard", "Goose", "Embryonated_eggs", "Sus_scrofa",
    "American_green-winged_teal", "American_wigeon", "Blue-winged_teal",
    "Green-winged_teal", "Greater_white-fronted_goose", "Spot-billed_duck",
    "Bufflehead", "Aquatic_bird", "Domestic_duck", "Golden_monkey",
}


def to_canonical(s: str) -> str:
    """Canonical form for matching."""
    s = s.replace(" ", "_").strip()
    if not s.startswith("A/"):
        return s
    parts = s.split("/")
    if parts[-1].endswith("(H3N2)"):
        subtype = "(H3N2)"
        parts[-1] = parts[-1].replace("(H3N2)", "")
    else:
        subtype = ""
    if len(parts) >= 3 and parts[1] in KNOWN_HOSTS:
        canonical = "A/" + "/".join(parts[2:]) + subtype
    else:
        canonical = "/".join(parts) + subtype
    return canonical


def parse_fasta_only_from_report() -> set[str]:
    """Parse FASTA ONLY strains from comparison report."""
    strains = set()
    if not REPORT_FILE.exists():
        print(f"Warning: {REPORT_FILE} not found. Run compare_fasta_metadata.py first.")
        return strains
    in_fasta_only = False
    with open(REPORT_FILE, encoding="utf-8") as f:
        for line in f:
            if "FASTA ONLY" in line:
                in_fasta_only = True
                continue
            if in_fasta_only and line.strip().startswith("A/"):
                strains.add(line.strip().strip())
    return strains


def load_bvbrc_data() -> list[dict]:
    """Load BV-BRC metadata from genome_metadata_H3N2.tsv or run p3-all-genomes."""
    if GENOME_META.exists():
        rows = []
        with open(GENOME_META, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            raw_rows = list(reader)
            # Normalize keys: genome.genome_id -> genome_id for compatibility
            for row in raw_rows:
                normalized = {}
                for k, v in row.items():
                    key = k.split(".")[-1] if "." in k else k
                    normalized[key] = v
                rows.append(normalized)
        print(f"Loaded {len(rows)} rows from {GENOME_META}")
        return rows

    # Try running p3-all-genomes
    print("genome_metadata_H3N2.tsv not found. Running p3-all-genomes...")
    try:
        result = subprocess.run(
            [
                "p3-all-genomes",
                '--eq=genome_name,Influenza A virus',
                "--eq=subtype,H3N2",
                "--attr", "genome_id,genome_name,strain,collection_date,isolation_country,geographic_group,host_common_name,segment",
            ],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode != 0:
            print(f"p3-all-genomes failed: {result.stderr}")
            return []
        # Parse TSV output
        rows = []
        lines = result.stdout.strip().split("\n")
        if not lines:
            return []
        header = lines[0].split("\t")
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= len(header):
                rows.append(dict(zip(header, parts)))
        print(f"Fetched {len(rows)} rows from p3-all-genomes")
        return rows
    except FileNotFoundError:
        print('p3-all-genomes not found. Run: p3-all-genomes --eq="genome_name,Influenza A virus" --eq=subtype,H3N2 ... > genome_metadata_H3N2.tsv')
        return []
    except subprocess.TimeoutExpired:
        print("p3-all-genomes timed out. Run manually and save to genome_metadata_H3N2.tsv")
        return []


def extract_strain_from_bvbrc(row: dict) -> str | None:
    """Extract strain ID from BV-BRC row (genome_name or strain column)."""
    gn = row.get("genome_name", row.get("genome name", "")).strip().strip('"')
    if gn:
        m = re.search(r"\(?(A/[^()]*(?:\(H3N2\))?)\)?\s*$", gn)
        if m:
            return m.group(1).strip()
        m2 = re.search(r"(A/[\w/_-]+(?:/\d{4})?)", gn)
        if m2:
            s = m2.group(1)
            return s + "(H3N2)" if "(H3N2)" not in s else s
    strain = row.get("strain", "").strip().strip('"')
    if strain and strain.startswith("A/"):
        return strain + "(H3N2)" if "(H3N2)" not in strain else strain
    return None


def main():
    fasta_only = parse_fasta_only_from_report()
    if not fasta_only:
        print("No FASTA-only strains found. Run compare_fasta_metadata.py first.")
        return 1

    print(f"FASTA-only strains: {len(fasta_only)}")
    fasta_canon = {to_canonical(s): s for s in fasta_only}

    bvbrc_rows = load_bvbrc_data()
    if not bvbrc_rows:
        print("No BV-BRC data in genome_metadata_H3N2.tsv (file empty or missing).")
        print('Re-run: p3-all-genomes --eq="genome_name,Influenza A virus" --eq=subtype,H3N2 --attr ... > genome_metadata_H3N2.tsv')
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        with open(ENRICHED_FILE, "w", encoding="utf-8", newline="") as f:
            f.write("strain\tcountry\tregion\thost\tnum_date\n")
        print(f"Wrote empty {ENRICHED_FILE}. Run tune_all_auspice_meta.py anyway (uses NEW_genome_w_clade + LOCATION_TO_COUNTRY).")
        return 0

    enriched = {}
    for row in bvbrc_rows:
        strain = extract_strain_from_bvbrc(row)
        if not strain:
            continue
        canon = to_canonical(strain)
        if canon not in fasta_canon:
            continue
        orig = fasta_canon[canon]
        if orig in enriched:
            continue
        country = row.get("isolation_country", row.get("isolation country", "")).strip().strip('"')
        region = row.get("geographic_group", row.get("geographic group", "")).strip().strip('"')
        host = row.get("host_common_name", row.get("host common name", "")).strip().strip('"')
        date = row.get("collection_date", row.get("collection date", "")).strip().strip('"')
        if country and not region:
            region = get_region_for_country(country) or ""
        # If geographic_group is a location (e.g. Hong_Kong), resolve to country
        geo = row.get("geographic_group", row.get("geographic group", "")).strip()
        if geo in LOCATION_TO_COUNTRY and not country:
            country, region = LOCATION_TO_COUNTRY[geo]
        num_date = None
        if date:
            m = re.match(r"(\d{4})", date)
            if m:
                num_date = float(m.group(1)) + 0.5
        enriched[orig] = {
            "strain": orig,
            "country": country,
            "region": region,
            "host": host,
            "num_date": num_date,
        }

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(ENRICHED_FILE, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["strain", "country", "region", "host", "num_date"], delimiter="\t")
        w.writeheader()
        w.writerows(enriched.values())

    print(f"Enriched {len(enriched)} strains -> {ENRICHED_FILE}")
    return 0


if __name__ == "__main__":
    exit(main())
