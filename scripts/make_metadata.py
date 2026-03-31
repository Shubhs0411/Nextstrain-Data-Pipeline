"""
Convert BV-BRC metadata to Nextstrain format.
Merge NEW_genome_w_clade.txt for country/region/host/clade.
Fill region from country when region is empty.
"""
import csv
import re
import sys
from pathlib import Path

# Add parent for location_to_country
sys.path.insert(0, str(Path(__file__).resolve().parent))
from location_to_country import get_region_for_country, LOCATION_TO_COUNTRY

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
OUTPUT_DIR = PROJECT_ROOT / "H3N2_output"
BVBRC_FILE = DATA_DIR / "BVBRC_genome.txt"
NEW_CLADE_FILE = PROJECT_ROOT / "NEW_genome_w_clade.txt"
OUTPUT_FILE = OUTPUT_DIR / "metadata.tsv"

# BV-BRC column names (from header)
BVBRC_COLS = [
    "Genome ID", "Genome Name", "Other Names", "NCBI Taxon ID", "Taxon Lineage IDs",
    "Taxon Lineage Names", "Superkingdom", "Kingdom", "Phylum", "Class", "Order",
    "Family", "Genus", "Species", "Genome Status", "Strain", "Serovar", "Biovar",
    "Pathovar", "MLST", "Segment", "Subtype", "H_type", "N_type", "H1 Clade Global",
    "H1 Clade US", "H5 Clade", "pH1N1-like", "Lineage", "Clade", "Subclade",
    "Other Typing", "Culture Collection", "Type Strain", "Reference", "Genome Quality",
    "Completion Date", "Publication", "Authors", "BioProject Accession",
    "BioSample Accession", "Assembly Accession", "SRA Accession", "GenBank Accessions",
    "Sequencing Center", "Sequencing Status", "Sequencing Platform", "Sequencing Depth",
    "Assembly Method", "Chromosome", "Plasmids", "Contigs", "Size", "GC Content",
    "Contig L50", "Contig N50", "TRNA", "RRNA", "Mat Peptide", "CDS",
    "Coarse Consistency", "Fine Consistency", "CheckM Contamination", "CheckM Completeness",
    "Genome Quality Flags", "Isolation Source", "Isolation Comments", "Collection Date",
    "Collection Year", "Season", "Isolation Country", "State/Province", "County",
    "City", "Geographic Group", "Geographic Location", "Other Environmental",
    "Host Name", "Host Common Name", "Host Sex", "Host Age", "Host Health", "Host Group",
    "Lab Host", "Passage", "Other Clinical", "Additional Metadata", "Comments",
    "Date Inserted", "Date Modified",
]


def strip(s: str) -> str:
    return s.strip().strip('"')


def normalize_date(s: str) -> str:
    """Normalize to YYYY-MM-DD or YYYY. Strip time."""
    s = strip(s)
    if not s:
        return ""
    # Remove time part
    s = re.sub(r"T\d{2}:\d{2}:\d{2}.*$", "", s)
    s = s.strip()
    if re.match(r"^\d{4}-\d{2}-\d{2}$", s):
        return s
    if re.match(r"^\d{4}-\d{2}$", s):
        return s + "-01"
    if re.match(r"^\d{4}$", s):
        return s + "-01-01"
    return s


def pubmed_url(pub: str) -> str:
    """Convert PMID to PubMed URL."""
    pub = strip(pub)
    if not pub:
        return ""
    m = re.search(r"\d+", pub)
    if m:
        return f"https://www.ncbi.nlm.nih.gov/pubmed/{m.group()}"
    return pub


def load_new_genome_w_clade() -> dict[str, dict]:
    """Load NEW_genome_w_clade.txt: strain_id -> {host, country, region, year, clade}."""
    result = {}
    if not NEW_CLADE_FILE.exists():
        return result
    with open(NEW_CLADE_FILE, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            strain_id = strip(parts[0])
            host = strip(parts[1]) if len(parts) > 1 else ""
            country_region = strip(parts[2]) if len(parts) > 2 else ""
            year = strip(parts[3]) if len(parts) > 3 else ""
            clade = strip(parts[4]) if len(parts) > 4 else ""  # may be "unassigned"

            # Normalize strain_id: (A/Human/...) -> A/Human/...
            clean = strain_id
            if clean.startswith("(") and clean.endswith(")"):
                clean = clean[1:-1]
            clean = clean.replace(" ", "_").replace("(H3N2)", "").strip()

            # Resolve country/region from NEW file column (can be country or location)
            country_val = country_region.replace(" ", "_")
            region_val = ""
            if country_val in LOCATION_TO_COUNTRY:
                country_val, region_val = LOCATION_TO_COUNTRY[country_val]
            elif country_val:
                region_val = get_region_for_country(country_val) or ""

            result[clean] = {
                "host": host,
                "country": country_val,
                "region": region_val,
                "year": year,
                "clade": clade,
            }
            # Also store with (H3N2) suffix for matching
            if "(H3N2)" not in clean:
                result[clean + "(H3N2)"] = result[clean].copy()
            # Also store without host for matching (A/Human/X -> A/X)
            if clean.startswith("A/Human/") and "/" in clean[8:]:
                alt = "A/" + clean[8:]
                if alt not in result:
                    result[alt] = result[clean].copy()
    return result


def main():
    new_clade = load_new_genome_w_clade()
    print(f"Loaded {len(new_clade)} strains from NEW_genome_w_clade.txt")

    with open(BVBRC_FILE, encoding="utf-8", errors="replace") as fin:
        header = fin.readline().strip().split("\t")
        col_idx = {h: i for i, h in enumerate(header)}

    def get(row: list[str], name: str) -> str:
        i = col_idx.get(name, -1)
        return strip(row[i]) if i >= 0 and i < len(row) else ""

    out_cols = [
        "strain", "virus", "accession", "date", "region", "country", "division", "city",
        "db", "segment", "host", "authors", "url", "title", "journal", "paper_url", "clade",
    ]

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(BVBRC_FILE, encoding="utf-8", errors="replace") as fin, \
         open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(out_cols)

        header = next(reader)
        col_idx = {h: i for i, h in enumerate(header)}

        for row in reader:
            if len(row) <= max(col_idx.get("Strain", 0), col_idx.get("Segment", 0)):
                continue

            strain = get(row, "Strain")
            if not strain or not strain.startswith("A/"):
                continue

            segment = get(row, "Segment")
            date = normalize_date(get(row, "Collection Date"))
            country = get(row, "Isolation Country")
            region = get(row, "Geographic Group")
            division = get(row, "State/Province")
            city = get(row, "City")
            host_common = get(row, "Host Common Name")
            host_name = get(row, "Host Name")
            host = host_common or host_name or ""
            authors = get(row, "Authors")
            pub = get(row, "Publication")
            accession = get(row, "GenBank Accessions")

            # Fill region from country when empty
            if country and not region:
                region = get_region_for_country(country) or ""

            # Merge NEW_genome_w_clade
            strain_clean = strain.replace(" ", "_")
            if strain_clean in new_clade:
                nc = new_clade[strain_clean]
                if nc.get("country"):
                    country = country or nc["country"]
                if nc.get("region"):
                    region = region or nc["region"]
                if nc.get("host"):
                    host = host or nc["host"]
                if nc.get("year") and not date:
                    date = nc["year"] + "-01-01" if len(nc["year"]) == 4 else nc["year"]

            clade = new_clade.get(strain_clean, {}).get("clade", "") if strain_clean in new_clade else ""

            writer.writerow([
                strain,
                "flu_h3n2",
                accession,
                date,
                region,
                country,
                division,
                city,
                "BV-BRC",
                segment,
                host,
                authors,
                pubmed_url(pub),
                "",
                "",
                pubmed_url(pub),
                clade,
            ])
            count += 1
            if count % 50000 == 0:
                print(f"  Processed {count} rows...")

    print(f"Wrote {count} rows to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
