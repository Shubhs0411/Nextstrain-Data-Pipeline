#!/usr/bin/env python3
"""
Build one Path C concatenated dataset (segments 1..8) and export a single Auspice JSON.

Outputs:
  - H3N2_DATA/H3N2_DATA/h3n2_pathc_concat.clean.fasta
  - H3N2_DATA/H3N2_DATA/metadata_pathc_concat.tsv
  - results/pathc_concat/{aligned.fasta,tree.nwk,tree-refined.nwk,traits_host.json}
  - auspice/h3n2_pathc_concat.json
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "H3N2_DATA" / "H3N2_DATA"
RESULTS_DIR = PROJECT_ROOT / "results" / "pathc_concat"
AUSPICE_DIR = PROJECT_ROOT / "auspice"
CONFIG = PROJECT_ROOT / "config" / "config.yaml"

SEGMENTS = list(range(1, 9))
# Merge demographics from HA first, then fill gaps from other segments.
MERGE_ORDER_DEMO = [4, 1, 2, 3, 5, 6, 7, 8]
# Clade labels: HA then NA, then remaining segments.
MERGE_ORDER_CLADE = [4, 6, 1, 2, 3, 5, 7, 8]


def run(cmd: list[str]) -> None:
    print(" ", " ".join(cmd))
    try:
        subprocess.run(cmd, cwd=str(PROJECT_ROOT), check=True)
    except FileNotFoundError as e:
        exe = cmd[0] if cmd else "?"
        raise RuntimeError(
            f"Could not run {exe!r} (not on PATH). On native Windows this often fails even if "
            f"Python works. Run the concat build in WSL with your Nextstrain env, e.g.:\n"
            f"  conda activate bvbrc\n"
            f"  cd /mnt/c/Users/sdeshmuk/Desktop/H3N2\n"
            f"  python scripts/build_pathc_concat.py"
        ) from e


def check_tree_build_tools() -> None:
    """augur align needs mafft; augur tree needs iqtree2 (or iqtree)."""
    missing: list[str] = []
    if not shutil.which("augur"):
        missing.append("augur")
    if not shutil.which("mafft"):
        missing.append("mafft")
    if not shutil.which("iqtree2") and not shutil.which("iqtree"):
        missing.append("iqtree2")
    if missing:
        print(
            "\nERROR: Missing on PATH: " + ", ".join(missing),
            file=sys.stderr,
        )
        print(
            "Build the concat tree in WSL (or Git Bash) where conda `bvbrc` provides augur/mafft/iqtree.",
            file=sys.stderr,
        )
        sys.exit(1)


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    if not path.exists():
        return seqs
    current = None
    chunks: list[str] = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None and current not in seqs:
                    seqs[current] = "".join(chunks).upper()
                current = line[1:].split()[0].strip()
                chunks = []
            else:
                chunks.append(line)
        if current is not None and current not in seqs:
            seqs[current] = "".join(chunks).upper()
    return seqs


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    if not path.exists():
        return out
    with open(path, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gid = (row.get("strain") or "").strip().strip('"')
            if gid:
                out[gid] = row
    return out


def parse_config() -> dict[str, str]:
    cfg: dict[str, str] = {}
    if not CONFIG.exists():
        return cfg
    with open(CONFIG, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.match(r"(\w+):\s*(.+)", line.strip())
            if m:
                cfg[m.group(1).strip()] = m.group(2).strip().strip('"')
    return cfg


def _strip_cell(v: str | None) -> str:
    return (v or "").strip().strip('"')


def _nonempty_clade(v: str) -> bool:
    s = _strip_cell(v)
    return bool(s) and s.lower() != "unassigned"


def collect_concat_fieldnames() -> list[str]:
    """Union of all segment metadata columns; HA (segment 4) column order first, then new keys."""
    seen: set[str] = set()
    out: list[str] = []
    for seg in MERGE_ORDER_DEMO:
        p = DATA_DIR / f"metadata_segment{seg}_pathc.tsv"
        with open(p, encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for h in reader.fieldnames or []:
                if h not in seen:
                    seen.add(h)
                    out.append(h)
    return out


def merge_row_for_strain(
    sname: str,
    by_name_by_seg: dict[int, dict[str, dict[str, str]]],
    fieldnames: list[str],
) -> dict[str, str]:
    """
    One row per isolate: fill every column from segment TSVs so nothing is dropped.
    Demographics: first non-empty across MERGE_ORDER_DEMO (HA first).
    clade/subclade: HA -> NA -> others; skip empty and 'unassigned'.
    """
    rows = {seg: by_name_by_seg[seg][sname] for seg in SEGMENTS}
    merged: dict[str, str] = {k: "" for k in fieldnames}

    for k in fieldnames:
        if k in ("strain", "segment", "clade", "subclade"):
            continue
        for seg in MERGE_ORDER_DEMO:
            v = _strip_cell(rows[seg].get(k))
            if v:
                merged[k] = v
                break

    merged["clade"] = ""
    merged["subclade"] = ""
    for seg in MERGE_ORDER_CLADE:
        r = rows[seg]
        c = _strip_cell(r.get("clade"))
        if not merged["clade"] and _nonempty_clade(c):
            merged["clade"] = c
        sc = _strip_cell(r.get("subclade"))
        if not merged["subclade"] and _nonempty_clade(sc):
            merged["subclade"] = sc

    return merged


def build_concat_inputs() -> tuple[Path, Path, int]:
    fasta_by_seg: dict[int, dict[str, str]] = {}
    meta_by_seg: dict[int, dict[str, dict[str, str]]] = {}
    by_name_by_seg: dict[int, dict[str, dict[str, str]]] = {}

    for seg in SEGMENTS:
        fasta_path = DATA_DIR / f"segment_{seg}_h3n2_pathc.clean.fasta"
        meta_path = DATA_DIR / f"metadata_segment{seg}_pathc.tsv"
        fseqs = read_fasta(fasta_path)
        mrows = read_metadata(meta_path)
        if not fseqs:
            raise RuntimeError(f"Missing/empty FASTA: {fasta_path}")
        if not mrows:
            raise RuntimeError(f"Missing/empty metadata: {meta_path}")
        fasta_by_seg[seg] = fseqs
        meta_by_seg[seg] = mrows
        by_name: dict[str, dict[str, str]] = {}
        for gid, row in mrows.items():
            sname = (row.get("strain_name") or "").strip().strip('"')
            if sname and sname not in by_name:
                by_name[sname] = row
        by_name_by_seg[seg] = by_name
        print(f"Segment {seg}: FASTA={len(fseqs)} metadata={len(mrows)}")

    # Path C genome_ids differ by segment. Join by strain_name, then map each segment
    # back to that segment's genome_id/sequence.
    common_names = set(by_name_by_seg[1].keys())
    for seg in SEGMENTS[1:]:
        common_names &= set(by_name_by_seg[seg].keys())
    common_names = sorted(common_names)
    if not common_names:
        raise RuntimeError("No common strain_name values across all 8 segment metadata files.")

    rows_by_name: list[tuple[str, dict[int, str]]] = []
    for sname in common_names:
        seg_to_gid: dict[int, str] = {}
        ok = True
        for seg in SEGMENTS:
            row = by_name_by_seg[seg].get(sname)
            gid = (row.get("strain") or "").strip().strip('"') if row else ""
            if not gid or gid not in fasta_by_seg[seg]:
                ok = False
                break
            seg_to_gid[seg] = gid
        if ok:
            rows_by_name.append((sname, seg_to_gid))
    if not rows_by_name:
        raise RuntimeError("No common strain_name values had all 8 segment sequences available.")

    concat_fasta = DATA_DIR / "h3n2_pathc_concat.clean.fasta"
    concat_meta = DATA_DIR / "metadata_pathc_concat.tsv"

    fieldnames = collect_concat_fieldnames()
    for req in ("strain", "strain_name", "segment", "clade", "subclade"):
        if req not in fieldnames:
            fieldnames.append(req)

    with open(concat_fasta, "w", encoding="utf-8") as ff, open(
        concat_meta, "w", encoding="utf-8", newline=""
    ) as fm:
        writer = csv.DictWriter(fm, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for sname, seg_to_gid in rows_by_name:
            pieces = [fasta_by_seg[seg][seg_to_gid[seg]] for seg in SEGMENTS]
            concat_seq = "".join(pieces)
            out_id = seg_to_gid[4]  # Keep HA genome_id as the node name in concat output.
            ff.write(f">{out_id}\n{concat_seq}\n")

            row = merge_row_for_strain(sname, by_name_by_seg, fieldnames)
            row["strain"] = out_id
            row["segment"] = "concat_1to8"
            writer.writerow(row)

    print(f"Wrote concatenated FASTA: {concat_fasta}")
    print(f"Wrote concatenated metadata: {concat_meta}")
    print(f"Common strains across segments: {len(rows_by_name)}")
    return concat_fasta, concat_meta, len(rows_by_name)


def build_concat_tree(concat_fasta: Path, concat_meta: Path, output_json: Path) -> None:
    cfg = parse_config()
    clock = cfg.get("clock_rate", "0.0008")
    coalescent = cfg.get("coalescent", "skyline")
    nthreads = cfg.get("align_nthreads", "auto")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    AUSPICE_DIR.mkdir(parents=True, exist_ok=True)

    aligned = RESULTS_DIR / "aligned.fasta"
    tree = RESULTS_DIR / "tree.nwk"
    tree_refined = RESULTS_DIR / "tree-refined.nwk"
    traits = RESULTS_DIR / "traits_host.json"

    print("1. Aligning concatenated sequences...")
    run(["augur", "align", "--sequences", str(concat_fasta), "--output", str(aligned), "--nthreads", str(nthreads)])

    print("2. Building tree...")
    run(["augur", "tree", "--alignment", str(aligned), "--output", str(tree)])

    print("3. Refining tree...")
    run(
        [
            "augur",
            "refine",
            "--tree",
            str(tree),
            "--alignment",
            str(aligned),
            "--metadata",
            str(concat_meta),
            "--output-tree",
            str(tree_refined),
            "--clock-rate",
            str(clock),
            "--coalescent",
            str(coalescent),
        ]
    )

    print("4. Inferring host traits...")
    run(
        [
            "augur",
            "traits",
            "--tree",
            str(tree_refined),
            "--metadata",
            str(concat_meta),
            "--output",
            str(traits),
            "--columns",
            "host",
        ]
    )

    print("5. Exporting Auspice JSON...")
    run(
        [
            "augur",
            "export",
            "v2",
            "--tree",
            str(tree_refined),
            "--metadata",
            str(concat_meta),
            "--node-data",
            str(traits),
            "--output",
            str(output_json),
        ]
    )

    geo_js = PROJECT_ROOT / "scripts" / "normalize_geo_resolutions.js"
    if geo_js.exists():
        print("6. Normalizing geo_resolutions (Auspice dir)...")
        run(["node", str(geo_js)])


def main() -> int:
    ap = argparse.ArgumentParser(description="Build one concatenated Path C Auspice JSON")
    ap.add_argument(
        "--output-json",
        default=str(AUSPICE_DIR / "h3n2_pathc_concat.json"),
        help="Output Auspice JSON path (default: auspice/h3n2_pathc_concat.json)",
    )
    args = ap.parse_args()

    output_json = Path(args.output_json)
    if not output_json.is_absolute():
        output_json = PROJECT_ROOT / output_json

    try:
        concat_fasta, concat_meta, n = build_concat_inputs()
        print(f"\nBuilding tree from {n} concatenated genomes...\n")
        check_tree_build_tools()
        build_concat_tree(concat_fasta, concat_meta, output_json)
    except subprocess.CalledProcessError as e:
        print(f"\nError: command failed ({e.returncode}): {' '.join(e.cmd)}", file=sys.stderr)
        return e.returncode or 1
    except RuntimeError as e:
        print(f"\nError: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        return 1

    print(f"\nDone. Output JSON: {output_json}")
    print("View with: npx auspice view --datasetDir auspice")
    return 0


if __name__ == "__main__":
    sys.exit(main())

