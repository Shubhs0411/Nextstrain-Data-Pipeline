"""
Clean FASTA headers for H3N2 segment files.
- Remove outer parentheses: (A/Human/...) -> A/Human/...
- Remove (H3N2) suffix: .../2013(H3N2) -> .../2013
- Normalize BV-BRC human naming: A/Human/... -> A/... (so FASTA matches BV-BRC Strain)
"""
import re
import sys
from pathlib import Path


def clean_header(header: str) -> str:
    """Clean a FASTA header line for Nextstrain/BV-BRC matching."""
    if not header.startswith('>'):
        return header
    h = header.strip().lstrip('>').strip()
    # Remove outer parentheses
    if h.startswith('(') and h.endswith(')'):
        h = h[1:-1]
    # Remove (H3N2) suffix
    h = re.sub(r'\(H3N2\)\s*$', '', h, flags=re.IGNORECASE).strip()
    # Normalize BV-BRC human naming: A/Human/... -> A/...
    if h.startswith('A/Human/'):
        h = 'A/' + h[8:]  # Remove "Human/"
    return '>' + h


def clean_fasta_file(input_path: Path, output_path: Path) -> int:
    """Clean headers in a FASTA file. Returns count of sequences processed."""
    count = 0
    with open(input_path, encoding='utf-8', errors='replace') as fin, \
         open(output_path, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('>'):
                fout.write(clean_header(line) + '\n')
                count += 1
            else:
                fout.write(line)
    return count


def main():
    if len(sys.argv) < 2:
        print("Usage: python clean_fasta_headers.py <input.fasta> [output.fasta]")
        sys.exit(1)
    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else input_path.with_stem(
        input_path.stem + '.clean'
    )
    if not input_path.exists():
        print(f"Error: {input_path} not found")
        sys.exit(1)
    count = clean_fasta_file(input_path, output_path)
    print(f"Cleaned {count} sequences -> {output_path}")


if __name__ == '__main__':
    main()
