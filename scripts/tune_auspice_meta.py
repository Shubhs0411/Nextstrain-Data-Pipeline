"""
Post-process Auspice JSON: inject metadata for tips, add colorings, filters, geo_resolutions.
- Lookup strain in metadata + NEW_genome_w_clade (bidirectional matching)
- Fallback: parse year/host from strain name; LOCATION_TO_COUNTRY inference for country/region
"""
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from location_to_country import infer_country_region_from_strain, get_region_for_country, LOCATION_TO_COUNTRY

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "H3N2_output"
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
AUSPICE_DIR = PROJECT_ROOT / "auspice"
NEW_CLADE_FILE = PROJECT_ROOT / "NEW_genome_w_clade.txt"

# Geo resolutions: country/region -> {latitude, longitude}
GEO_RESOLUTIONS = {
    "country": {
        "USA": {"latitude": 38.0, "longitude": -97.0},
        "China": {"latitude": 35.0, "longitude": 105.0},
        "Japan": {"latitude": 36.0, "longitude": 138.0},
        "Australia": {"latitude": -25.0, "longitude": 133.0},
        "United Kingdom": {"latitude": 54.0, "longitude": -2.0},
        "Germany": {"latitude": 51.0, "longitude": 10.0},
        "France": {"latitude": 46.0, "longitude": 2.0},
        "Italy": {"latitude": 43.0, "longitude": 12.0},
        "Spain": {"latitude": 40.0, "longitude": -4.0},
        "Canada": {"latitude": 56.0, "longitude": -106.0},
        "Mexico": {"latitude": 23.0, "longitude": -102.0},
        "Brazil": {"latitude": -10.0, "longitude": -55.0},
        "India": {"latitude": 20.0, "longitude": 77.0},
        "South Korea": {"latitude": 36.0, "longitude": 128.0},
        "Thailand": {"latitude": 15.0, "longitude": 101.0},
        "Vietnam": {"latitude": 16.0, "longitude": 108.0},
        "Singapore": {"latitude": 1.3, "longitude": 103.8},
        "Denmark": {"latitude": 56.0, "longitude": 10.0},
        "Netherlands": {"latitude": 52.0, "longitude": 5.0},
        "Sweden": {"latitude": 62.0, "longitude": 15.0},
        "Norway": {"latitude": 62.0, "longitude": 10.0},
        "New Zealand": {"latitude": -41.0, "longitude": 174.0},
        "Chile": {"latitude": -35.0, "longitude": -71.0},
        "Argentina": {"latitude": -34.0, "longitude": -64.0},
        "Russia": {"latitude": 60.0, "longitude": 100.0},
        "Mongolia": {"latitude": 46.0, "longitude": 105.0},
        "Iran": {"latitude": 32.0, "longitude": 53.0},
        "Bangladesh": {"latitude": 24.0, "longitude": 90.0},
        "Indonesia": {"latitude": -5.0, "longitude": 120.0},
        "Malaysia": {"latitude": 4.0, "longitude": 102.0},
        "Philippines": {"latitude": 13.0, "longitude": 122.0},
        "Cambodia": {"latitude": 13.0, "longitude": 105.0},
        "Egypt": {"latitude": 27.0, "longitude": 30.0},
        "Uganda": {"latitude": 1.0, "longitude": 32.0},
        "Kenya": {"latitude": 0.0, "longitude": 38.0},
        "South Africa": {"latitude": -29.0, "longitude": 24.0},
        "Greece": {"latitude": 39.0, "longitude": 22.0},
        "Turkey": {"latitude": 39.0, "longitude": 35.0},
        "Nicaragua": {"latitude": 13.0, "longitude": -85.0},
        "Panama": {"latitude": 9.0, "longitude": -80.0},
        "Colombia": {"latitude": 4.0, "longitude": -72.0},
        "Peru": {"latitude": -10.0, "longitude": -76.0},
        "Ecuador": {"latitude": -2.0, "longitude": -78.0},
        "Venezuela": {"latitude": 6.0, "longitude": -66.0},
        "Cuba": {"latitude": 22.0, "longitude": -80.0},
        "Puerto Rico": {"latitude": 18.0, "longitude": -66.0},
        "Dominican Republic": {"latitude": 19.0, "longitude": -70.0},
        "Haiti": {"latitude": 19.0, "longitude": -72.0},
        "El Salvador": {"latitude": 13.8, "longitude": -88.9},
        "Guatemala": {"latitude": 15.5, "longitude": -90.3},
        "Honduras": {"latitude": 15.2, "longitude": -86.2},
        "Costa Rica": {"latitude": 10.0, "longitude": -84.0},
        "Jamaica": {"latitude": 18.2, "longitude": -77.3},
    },
    "region": {
        "North America": {"latitude": 38.0, "longitude": -97.0},
        "South America": {"latitude": -14.0, "longitude": -59.0},
        "Europe": {"latitude": 50.0, "longitude": 10.0},
        "Asia": {"latitude": 34.0, "longitude": 105.0},
        "Africa": {"latitude": 4.0, "longitude": 22.0},
        "Oceania": {"latitude": -25.0, "longitude": 133.0},
    },
}


def parse_year_from_strain(strain: str) -> float | None:
    """Extract year from strain name (e.g. .../2013 -> 2013.5)."""
    m = re.search(r"/(\d{4})(?:\(H3N2\))?\s*$", strain)
    if m:
        return float(m.group(1)) + 0.5
    return None


def parse_host_from_strain(strain: str) -> str:
    """Extract host from strain name (e.g. A/Human/... -> Human)."""
    parts = strain.split("/")
    if len(parts) >= 3:
        host = parts[1]
        known = {"Human", "Swine", "Duck", "Chicken", "Canine", "Turkey", "Environment",
                 "Feline", "Mallard", "Goose", "Sus_scrofa", "Domestic_duck"}
        if host in known or host.replace("_", " ").title() in known:
            return host.replace("_", " ")
    return ""


def load_metadata_lookup(metadata_path: Path) -> dict[str, dict]:
    """Load metadata TSV into strain -> core fields plus any extra columns we know how to use."""
    lookup = {}
    if not metadata_path.exists():
        return lookup
    with open(metadata_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            strain = row.get("strain", "").strip().strip('"')
            if not strain:
                continue
            date = row.get("date", "").strip()
            num_date = None
            if date:
                m = re.match(r"(\d{4})", date)
                if m:
                    num_date = float(m.group(1)) + 0.5
            info = {
                "region": row.get("region", "").strip().strip('"'),
                "country": row.get("country", "").strip().strip('"'),
                "host": row.get("host", "").strip().strip('"'),
                "num_date": num_date,
                "clade": row.get("clade", "").strip().strip('"'),
            }
            # Optional richer metadata columns (mostly for Path C from BVBRC)
            for key in [
                "accession",
                "authors",
                "paper_url",
                "journal",
                "division",
                "city",
                "isolation_source",
                "subclade",
                "lineage",
            ]:
                if key in row:
                    info[key] = (row.get(key) or "").strip().strip('"')
            lookup[strain] = info
            # Bidirectional: add A/... form for A/Human/...
            if strain.startswith("A/Human/") and "/" in strain[8:]:
                alt = "A/" + strain[8:]
                if alt not in lookup:
                    lookup[alt] = lookup[strain].copy()
    return lookup


def load_bvbrc_enriched_lookup() -> dict[str, dict]:
    """Load bvbrc_enriched.tsv (from fetch_bvbrc_metadata.py) for FASTA-only strains."""
    path = OUTPUT_DIR / "bvbrc_enriched.tsv"
    lookup = {}
    if not path.exists():
        return lookup
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            strain = row.get("strain", "").strip()
            if not strain:
                continue
            num_date = row.get("num_date", "")
            if num_date:
                try:
                    num_date = float(num_date)
                except ValueError:
                    num_date = None
            lookup[strain] = {
                "country": row.get("country", "").strip(),
                "region": row.get("region", "").strip(),
                "host": row.get("host", "").strip(),
                "num_date": num_date,
                "clade": "",
            }
            if strain.startswith("A/Human/") and "/" in strain[8:]:
                alt = "A/" + strain[8:]
                if alt not in lookup:
                    lookup[alt] = lookup[strain].copy()
    return lookup


def load_new_clade_lookup() -> dict[str, dict]:
    """Load NEW_genome_w_clade.txt into strain -> {country, region, host, num_date, clade}."""
    lookup = {}
    if not NEW_CLADE_FILE.exists():
        return lookup
    with open(NEW_CLADE_FILE, encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            strain_id = parts[0].strip().strip('"')
            clean = strain_id
            if clean.startswith("(") and clean.endswith(")"):
                clean = clean[1:-1]
            clean = clean.replace(" ", "_").replace("(H3N2)", "").strip()
            host = parts[1].strip() if len(parts) > 1 else ""
            country_region = parts[2].strip().replace(" ", "_") if len(parts) > 2 else ""
            year = parts[3].strip() if len(parts) > 3 else ""
            clade = parts[4].strip() if len(parts) > 4 else ""

            country = country_region
            region = ""
            if country_region in LOCATION_TO_COUNTRY:
                country, region = LOCATION_TO_COUNTRY[country_region]
            elif country_region:
                region = get_region_for_country(country_region) or ""

            num_date = float(year) + 0.5 if year and year.isdigit() else None
            lookup[clean] = {"country": country, "region": region, "host": host, "num_date": num_date, "clade": clade}
            if clean.startswith("A/Human/") and "/" in clean[8:]:
                alt = "A/" + clean[8:]
                if alt not in lookup:
                    lookup[alt] = lookup[clean].copy()
    return lookup


def inject_node_attrs(node: dict, meta_lookup: dict, new_clade: dict, bvbrc_enriched: dict, all_countries: set, all_regions: set, all_hosts: set | None = None):
    """Recursively inject node_attrs for tips."""
    if all_hosts is None:
        all_hosts = set()
    name = node.get("name", "")
    if "children" in node:
        for child in node["children"]:
            inject_node_attrs(child, meta_lookup, new_clade, bvbrc_enriched, all_countries, all_regions, all_hosts)
        return

    # Leaf node
    if not name or name.startswith("NODE_"):
        return

    strain = name.strip()
    strain_alt = strain
    if strain.startswith("A/Human/") and "/" in strain[8:]:
        strain_alt = "A/" + strain[8:]
    elif "/" in strain and not strain.startswith("A/Human/"):
        strain_alt = "A/Human/" + strain[2:] if strain.startswith("A/") else strain

    attrs = node.setdefault("node_attrs", {})
    attrs["strain"] = {"value": strain}

    # Always clear any existing clade / subclade first so reruns don't keep
    # stale or 'unassigned' values.
    if "clade" in attrs:
        attrs.pop("clade", None)
    if "subclade" in attrs:
        attrs.pop("subclade", None)

    # Lookup: metadata first, then bvbrc_enriched, then NEW_genome_w_clade
    info = (
        meta_lookup.get(strain) or meta_lookup.get(strain_alt)
        or bvbrc_enriched.get(strain) or bvbrc_enriched.get(strain_alt)
        or new_clade.get(strain) or new_clade.get(strain_alt)
    )

    if info:
        if info.get("country"):
            attrs["country"] = {"value": info["country"]}
            all_countries.add(info["country"])
        if info.get("region"):
            attrs["region"] = {"value": info["region"]}
            all_regions.add(info["region"])
        if info.get("host"):
            attrs["host"] = {"value": info["host"]}
            all_hosts.add(info["host"])
        if info.get("num_date") is not None:
            attrs["num_date"] = {"value": info["num_date"]}
        clade = (info.get("clade") or "").strip().strip('"')
        if clade and clade.lower() != "unassigned":
            attrs["clade"] = {"value": clade}
        subclade = (info.get("subclade") or "").strip().strip('"')
        if subclade and subclade.lower() != "unassigned":
            attrs["subclade"] = {"value": subclade}
        if info.get("lineage"):
            attrs["lineage"] = {"value": info["lineage"]}
        # Extra optional fields (mainly from Path C / BV-BRC)
        if info.get("accession"):
            attrs["accession"] = {"value": info["accession"]}
        if info.get("authors"):
            attrs["authors"] = {"value": info["authors"]}
        if info.get("paper_url"):
            attrs["paper_url"] = {"value": info["paper_url"]}
        if info.get("journal"):
            attrs["journal"] = {"value": info["journal"]}
        if info.get("division"):
            attrs["division"] = {"value": info["division"]}
        if info.get("city"):
            attrs["city"] = {"value": info["city"]}
        if info.get("isolation_source"):
            attrs["isolation_source"] = {"value": info["isolation_source"]}
    else:
        # Fallback: parse from strain name
        num_date = parse_year_from_strain(strain)
        if num_date is not None:
            attrs["num_date"] = {"value": num_date}
        host = parse_host_from_strain(strain)
        if host:
            attrs["host"] = {"value": host}
            all_hosts.add(host)
        # LOCATION_TO_COUNTRY inference
        country, region = infer_country_region_from_strain(strain)
        if country:
            attrs["country"] = {"value": country}
            all_countries.add(country)
        if region:
            attrs["region"] = {"value": region}
            all_regions.add(region)


def ensure_geo_resolutions(meta: dict, all_countries: set, all_regions: set):
    """Add geo_resolutions for all countries/regions present."""
    gr = meta.get("geo_resolutions")
    # If geo_resolutions is already in Nextstrain v2 array form (from
    # normalize_geo_resolutions.js), don't try to treat it as a dict.
    if isinstance(gr, list):
        return
    # Otherwise, build / update legacy dict form and let the normalizer turn
    # it into an array later.
    gr = meta.setdefault("geo_resolutions", {})
    if "country" not in gr:
        gr["country"] = {"demes": {}}
    if "region" not in gr:
        gr["region"] = {"demes": {}}
    for c in all_countries:
        if c not in gr["country"]["demes"] and c in GEO_RESOLUTIONS["country"]:
            gr["country"]["demes"][c] = GEO_RESOLUTIONS["country"][c]
        elif c not in gr["country"]["demes"]:
            gr["country"]["demes"][c] = {"latitude": 0, "longitude": 0}
    for r in all_regions:
        if r not in gr["region"]["demes"] and r in GEO_RESOLUTIONS["region"]:
            gr["region"]["demes"][r] = GEO_RESOLUTIONS["region"][r]
        elif r not in gr["region"]["demes"]:
            gr["region"]["demes"][r] = {"latitude": 0, "longitude": 0}


# Nextstrain-style categorical color palette (hex); cycle for many values
CATEGORICAL_PALETTE = [
    "#447CCD", "#3F45C8", "#3E5DD0", "#4068CF", "#4B8EC1", "#549DB3", "#5EA8A2", "#6AB090",
    "#78B67E", "#8FBC66", "#A8BD54", "#B8BC4A", "#C8B944", "#D5B13F", "#E59638", "#E67531",
    "#E2592C", "#E04B29", "#DC2F24", "#4530BB", "#3F51CC", "#447DCC", "#58A3AB", "#64AC99",
    "#7FB975", "#97BD5F", "#B0BD4E", "#C0BA47", "#CFB541", "#E29E39", "#E68033", "#E4672E",
]


def _make_scale(values: set, palette: list[str] | None = None) -> list[list[str]]:
    """Build Zika-style scale: [[value, hexColor], ...] for sorted values."""
    palette = palette or CATEGORICAL_PALETTE
    out = []
    for i, v in enumerate(sorted(v for v in values if v)):
        out.append([v, palette[i % len(palette)]])
    return out


def enrich_meta_block(
    meta: dict,
    all_countries: set | None = None,
    all_regions: set | None = None,
    all_hosts: set | None = None,
):
    """Add colorings (with scale), filters, panels, display_defaults, and Nextstrain-style meta fields."""
    meta.setdefault("title", "H3N2 Influenza Phylogeny")
    meta.setdefault("updated", "2026-03-01")
    meta.setdefault("display_defaults", {})["tip_label"] = "strain"
    meta.setdefault("display_defaults", {})["map_triplicate"] = True
    # Nextstrain-style provenance (like Zika)
    meta.setdefault("build_url", "https://github.com/nextstrain/flu")
    meta.setdefault("data_provenance", [{"name": "BV-BRC", "url": "https://www.bv-brc.org/"}])
    meta.setdefault("maintainers", [{"name": "H3N2 build", "url": "https://nextstrain.org"}])
    meta.setdefault(
        "description",
        "H3N2 influenza segment phylogeny. Data from BV-BRC (Path A: genome_id-based). "
        "We acknowledge data generators and BV-BRC for open sharing of genetic data.",
    )
    colorings = meta.get("colorings", [])
    keys = {c["key"] for c in colorings if isinstance(c, dict)}
    for key in ["num_date", "country", "region", "host", "clade", "subclade"]:
        if key not in keys:
            c = {"key": key, "title": key.replace("_", " "), "type": "categorical" if key != "num_date" else "continuous"}
            colorings.append(c)
    # Add explicit color scale for categorical colorings (Zika-style)
    if all_countries:
        scale = _make_scale(all_countries)
        if scale:
            for c in colorings:
                if isinstance(c, dict) and c.get("key") == "country":
                    c["scale"] = scale
                    break
    if all_regions:
        scale = _make_scale(all_regions)
        if scale:
            for c in colorings:
                if isinstance(c, dict) and c.get("key") == "region":
                    c["scale"] = scale
                    break
    if all_hosts:
        scale = _make_scale(all_hosts)
        if scale:
            for c in colorings:
                if isinstance(c, dict) and c.get("key") == "host":
                    c["scale"] = scale
                    break
    meta["colorings"] = colorings
    meta.setdefault("filters", ["country", "region", "host"])
    if "panels" not in meta or not meta["panels"]:
        meta["panels"] = ["tree", "map"]
    elif "map" not in meta["panels"]:
        meta["panels"] = list(meta["panels"]) + ["map"]


def process_json(
    json_path: Path,
    segment: int,
    meta_suffix: str = "",
    meta_override: Path | None = None,
) -> bool:
    """Process a single Auspice JSON file. meta_suffix e.g. '_pathb' for Path B small build."""
    if meta_override is not None:
        meta_path = meta_override
        if not meta_path.exists():
            print(f"Metadata not found: {meta_path}")
            return False
    elif meta_suffix:
        meta_path = DATA_DIR / f"metadata_segment{segment}{meta_suffix}.tsv"
    else:
        meta_path = DATA_DIR / f"metadata_segment{segment}_h3n2.tsv"
    if meta_override is None and not meta_path.exists():
        meta_path = OUTPUT_DIR / f"metadata_segment{segment}.tsv"
    meta_lookup = load_metadata_lookup(meta_path)
    new_clade = load_new_clade_lookup()
    bvbrc_enriched = load_bvbrc_enriched_lookup()

    with open(json_path, encoding="utf-8") as f:
        data = json.load(f)

    all_countries = set()
    all_regions = set()
    all_hosts = set()

    if "tree" in data:
        inject_node_attrs(data["tree"], meta_lookup, new_clade, bvbrc_enriched, all_countries, all_regions, all_hosts)

    meta = data.setdefault("meta", {})
    enrich_meta_block(meta, all_countries=all_countries, all_regions=all_regions, all_hosts=all_hosts)
    ensure_geo_resolutions(meta, all_countries, all_regions)

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(f"Processed {json_path.name} ({len(all_countries)} countries, {len(all_regions)} regions)")
    return True


def main():
    if len(sys.argv) >= 2:
        arg = sys.argv[1]
        if arg.endswith(".json"):
            json_path = Path(arg)
            if not json_path.is_absolute():
                json_path = AUSPICE_DIR / json_path.name if (AUSPICE_DIR / json_path.name).exists() else json_path
            m = re.search(r"segment(\d+)", json_path.stem, re.I)
            seg = int(m.group(1)) if m else 1
            if "_pathc" in json_path.stem:
                meta_suffix = "_pathc"
            elif "_pathb" in json_path.stem:
                meta_suffix = "_pathb"
            else:
                meta_suffix = ""
            if json_path.exists():
                stem_l = json_path.stem.lower()
                if "pathc_concat" in stem_l:
                    process_json(
                        json_path,
                        seg,
                        meta_override=DATA_DIR / "metadata_pathc_concat.tsv",
                    )
                else:
                    process_json(json_path, seg, meta_suffix=meta_suffix)
            else:
                print(f"File not found: {json_path}")
        else:
            print("Usage: python tune_auspice_meta.py [h3n2_segmentN.json]")
    else:
        # Process all
        for seg in range(1, 9):
            jp = AUSPICE_DIR / f"h3n2_segment{seg}.json"
            if jp.exists():
                process_json(jp, seg)
            else:
                print(f"  Skip (not found): {jp}")


if __name__ == "__main__":
    main()
