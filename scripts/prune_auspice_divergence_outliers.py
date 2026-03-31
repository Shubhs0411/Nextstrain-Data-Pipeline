#!/usr/bin/env python3
"""Remove tip nodes from an Auspice v2 JSON where cumulative node_attrs.div exceeds a threshold."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path
from typing import Any


def get_div(node: dict) -> float:
    raw = (node.get("node_attrs") or {}).get("div")
    if isinstance(raw, dict):
        raw = raw.get("value", 0)
    try:
        return float(raw) if raw is not None else 0.0
    except (TypeError, ValueError):
        return 0.0


def is_leaf(node: dict) -> bool:
    return not node.get("children")


def prune_tree(node: dict, threshold: float) -> dict | None:
    if is_leaf(node):
        return None if get_div(node) > threshold else node

    children = node.get("children") or []
    kept: list[dict] = []
    for ch in children:
        pr = prune_tree(ch, threshold)
        if pr is not None:
            kept.append(pr)

    if not kept:
        return None
    if len(kept) == 1:
        return kept[0]

    out = node
    out["children"] = kept
    return out


def count_tips(node: dict) -> int:
    if is_leaf(node):
        return 1
    return sum(count_tips(c) for c in node.get("children") or [])


def collect_tip_divs(node: dict) -> list[float]:
    if is_leaf(node):
        return [get_div(node)]
    out: list[float] = []
    for c in node.get("children") or []:
        out.extend(collect_tip_divs(c))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Prune Auspice tips with div > threshold")
    ap.add_argument("json_path", type=Path, help="Input Auspice JSON")
    ap.add_argument("--threshold", type=float, default=0.6, help="Max cumulative div to keep (default 0.6)")
    ap.add_argument("--backup-suffix", type=str, default=".before_div_prune.json", help="Backup copy suffix")
    ap.add_argument("--no-backup", action="store_true")
    args = ap.parse_args()

    path = args.json_path
    if not path.exists():
        print(f"Not found: {path}", file=sys.stderr)
        return 1

    with path.open(encoding="utf-8") as f:
        data: dict[str, Any] = json.load(f)

    tree = data.get("tree")
    if not isinstance(tree, dict):
        print("Invalid JSON: missing tree", file=sys.stderr)
        return 1

    before = count_tips(tree)
    divs = collect_tip_divs(tree)
    outliers = sum(1 for d in divs if d > args.threshold)

    if not args.no_backup:
        backup = path.with_name(path.stem + args.backup_suffix)
        shutil.copy2(path, backup)
        print(f"Backup: {backup}")

    new_tree = prune_tree(tree, args.threshold)
    if new_tree is None:
        print("Error: entire tree pruned", file=sys.stderr)
        return 1

    after = count_tips(new_tree)
    data["tree"] = new_tree

    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    print(f"Tips before/after: {before} -> {after} (removed {before - after}, tips with div>{args.threshold} was {outliers})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
