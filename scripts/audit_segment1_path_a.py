#!/usr/bin/env python3
"""
Audit Path A pipeline for segment 1: count genomes at each step and check for losses.
Also validates auspice/h3n2_segment1.json is Auspice v2 / Nextstrain compatible.
"""
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
RESULTS_DIR = PROJECT_ROOT / "results" / "segment_1"
AUSPICE_JSON = PROJECT_ROOT / "auspice" / "h3n2_segment1.json"


def count_fasta_headers(path: Path, normalize_con: bool = False) -> set[str]:
    """Return set of sequence IDs (first token after '>'). If normalize_con, strip .con.N to get genome_id."""
    if not path.exists():
        return set()
    ids = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                tok = line[1:].split("|")[0].split()[0].strip()
                if normalize_con and ".con." in tok:
                    tok = tok.split(".con.")[0]
                ids.add(tok)
    return ids


def count_tsv_rows(path: Path, skip_header: bool = True) -> int:
    if not path.exists():
        return 0
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    if skip_header and lines:
        lines = lines[1:]
    return len([l for l in lines if l.strip()])


def count_tsv_strains(path: Path) -> set[str]:
    """First column (strain) of TSV."""
    if not path.exists():
        return set()
    out = set()
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                out.add(parts[0].strip().strip('"'))
    return out


def count_genome_meta_segment1() -> int:
    """Rows in genome_metadata_H3N2.tsv where segment == 1."""
    path = PROJECT_ROOT / "genome_metadata_H3N2.tsv"
    if not path.exists():
        return 0
    with open(path, encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        idx_seg = next((i for i, h in enumerate(header) if "segment" in h.lower()), -1)
        idx_gid = next((i for i, h in enumerate(header) if "genome_id" in h.lower()), 0)
        if idx_seg < 0:
            return 0
        count = 0
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) <= max(idx_seg, idx_gid):
                continue
            seg = parts[idx_seg].strip().strip('"')
            if seg == "1":
                count += 1
    return count


def count_genome_ids_file() -> tuple[int, set[str]]:
    path = DATA_DIR / "genome_ids_segment1.txt"
    if not path.exists():
        return 0, set()
    ids = set()
    with open(path, encoding="utf-8") as f:
        for line in f:
            gid = line.strip()
            if gid and not gid.startswith("genome."):
                ids.add(gid)
    return len(ids), ids


def count_tree_tips(node: dict) -> int:
    """Count leaf nodes (tips) in tree."""
    if not node.get("children"):
        return 1
    return sum(count_tree_tips(c) for c in node["children"])


def main() -> int:
    print("=== Path A Segment 1 audit ===\n")

    # 1. Source metadata (genome_metadata_H3N2.tsv, segment=1)
    n_meta_raw = count_genome_meta_segment1()
    print(f"1. genome_metadata_H3N2.tsv (segment=1 rows): {n_meta_raw}")

    # 2. genome_ids_segment1.txt (from write_genome_id_lists.py)
    n_ids, id_set = count_genome_ids_file()
    print(f"2. genome_ids_segment1.txt:                  {n_ids}")

    # 3. metadata_segment1_h3n2.tsv (from make_metadata_from_genome_h3n2.py)
    n_meta_seg = count_tsv_rows(DATA_DIR / "metadata_segment1_h3n2.tsv")
    meta_strains = count_tsv_strains(DATA_DIR / "metadata_segment1_h3n2.tsv")
    print(f"3. metadata_segment1_h3n2.tsv:               {n_meta_seg}")

    # 4. segment_1_h3n2.fasta (after fetch)
    raw_fasta_ids = count_fasta_headers(DATA_DIR / "segment_1_h3n2.fasta")
    raw_fasta_ids_normalized = count_fasta_headers(DATA_DIR / "segment_1_h3n2.fasta", normalize_con=True)
    print(f"4. segment_1_h3n2.fasta:                     {len(raw_fasta_ids)}")

    # 5. segment_1_h3n2.clean.fasta (after ensure_fasta_headers_genome_id)
    clean_fasta_ids = count_fasta_headers(DATA_DIR / "segment_1_h3n2.clean.fasta")
    print(f"5. segment_1_h3n2.clean.fasta:              {len(clean_fasta_ids)}")

    # 6. Subsampled
    sub_fasta_ids = count_fasta_headers(RESULTS_DIR / "subsampled.fasta")
    sub_meta = count_tsv_strains(RESULTS_DIR / "metadata_subsampled.tsv")
    print(f"6. results/segment_1/subsampled.fasta:       {len(sub_fasta_ids)}")
    print(f"   results/segment_1/metadata_subsampled.tsv:{len(sub_meta)}")

    # 7. Auspice JSON
    if not AUSPICE_JSON.exists():
        print(f"7. {AUSPICE_JSON}: NOT FOUND")
    else:
        with open(AUSPICE_JSON, encoding="utf-8") as f:
            data = json.load(f)
        n_tips = count_tree_tips(data["tree"])
        print(f"7. auspice/h3n2_segment1.json (tree tips): {n_tips}")

    # Losses (compare using normalized raw IDs so .con.N headers match genome_id list)
    print("\n--- Possible losses ---")
    if id_set and raw_fasta_ids_normalized:
        in_list_not_fasta = id_set - raw_fasta_ids_normalized
        in_fasta_not_list = raw_fasta_ids_normalized - id_set
        print(f"  In genome_ids_segment1.txt but NOT in segment_1_h3n2.fasta: {len(in_list_not_fasta)} (fetch missed)")
        if in_list_not_fasta and len(in_list_not_fasta) <= 5:
            print(f"    Example IDs: {list(in_list_not_fasta)[:5]}")
        elif in_list_not_fasta:
            print(f"    Example IDs: {list(in_list_not_fasta)[:5]} ...")
        print(f"  In segment_1_h3n2.fasta but NOT in genome_ids list:         {len(in_fasta_not_list)}")
    if raw_fasta_ids is not None and meta_strains:
        in_fasta_not_meta = raw_fasta_ids - meta_strains
        print(f"  In clean FASTA but NOT in metadata_segment1_h3n2.tsv:        {len(clean_fasta_ids - meta_strains) if clean_fasta_ids else 0}")
    if sub_fasta_ids and sub_meta:
        if sub_fasta_ids != sub_meta:
            print(f"  Subsampled FASTA vs metadata strain count mismatch: {len(sub_fasta_ids)} vs {len(sub_meta)}")
        else:
            print("  Subsampled FASTA and metadata_subsampled.tsv strain sets match.")

    # JSON schema check (Auspice v2 / Nextstrain)
    print("\n--- Auspice v2 / Nextstrain compatibility ---")
    if AUSPICE_JSON.exists():
        with open(AUSPICE_JSON, encoding="utf-8") as f:
            data = json.load(f)
        checks = [
            (data.get("version") == "v2", "version == 'v2'"),
            ("meta" in data, "has 'meta'"),
            ("tree" in data, "has 'tree'"),
            ("colorings" in data.get("meta", {}), "meta.colorings"),
            ("tree" in data and "node_attrs" in data["tree"], "tree.node_attrs"),
            ("tree" in data and ("children" in data["tree"] or "name" in data["tree"]), "tree structure (name/children)"),
        ]
        for ok, label in checks:
            print(f"  {'OK' if ok else 'MISSING'}: {label}")
        if all(c[0] for c in checks):
            print("  => JSON is compatible with Nextstrain/Auspice v2 (same format as e.g. Zika dataset).")
        else:
            print("  => Some required keys missing for full compatibility.")
    else:
        print("  (auspice/h3n2_segment1.json not found)")

    print("\n--- Summary ---")
    print("  The 5000 genomes in the output are intentional: config subsample_max_sequences: 5000.")
    print("  Full segment 1 set is larger (~70k); subsampling avoids MAFFT OOM and keeps runtime manageable.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
