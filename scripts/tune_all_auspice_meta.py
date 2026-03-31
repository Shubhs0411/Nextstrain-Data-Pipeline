"""
Run tune_auspice_meta.py on all segment JSON files in auspice/.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tune_auspice_meta import process_json

PROJECT_ROOT = Path(__file__).resolve().parent.parent
AUSPICE_DIR = PROJECT_ROOT / "auspice"


def main():
    count = 0
    for seg in range(1, 9):
        jp = AUSPICE_DIR / f"h3n2_segment{seg}.json"
        if jp.exists():
            process_json(jp, seg)
            count += 1
        else:
            print(f"  Skip (not found): {jp}")
    if count == 0:
        print("No Auspice JSON files found. Run Snakemake first to generate them.")
    else:
        print(f"Processed {count} JSON file(s).")


if __name__ == "__main__":
    main()
