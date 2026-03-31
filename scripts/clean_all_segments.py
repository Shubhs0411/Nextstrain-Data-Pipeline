"""
Clean FASTA headers for all H3N2 segment files (1-8).
Output: segment_1.clean.fasta ... segment_8.clean.fasta
"""
import sys
from pathlib import Path

# Add parent to path for clean_fasta_headers
sys.path.insert(0, str(Path(__file__).resolve().parent))
from clean_fasta_headers import clean_fasta_file

DATA_DIR = Path(__file__).resolve().parent.parent / "H3N2_DATA" / "H3N2_DATA"
OUTPUT_DIR = Path(__file__).resolve().parent.parent / "H3N2_output"
SEGMENTS = list(range(1, 9))


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for seg in SEGMENTS:
        input_path = DATA_DIR / f"segment_{seg}.fasta"
        output_path = OUTPUT_DIR / f"segment_{seg}.clean.fasta"
        if not input_path.exists():
            print(f"  WARNING: {input_path} not found, skipping")
            continue
        count = clean_fasta_file(input_path, output_path)
        print(f"Segment {seg}: cleaned {count} sequences -> {output_path}")


if __name__ == '__main__':
    main()
