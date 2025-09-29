#!/usr/bin/env python3
"""
Generate a FASTA file containing all single point mutants (19 per site, unless omitted) for
specified positions in a sequence.

Usage example:
  python point_mutant_fasta_cli.py \
      --sequence KKEEAMELVRKAASLLSPSKAEKVIAEVEAALEEGEEEAKKKLEEWARKAGAMKYSRAKKLLRKAAETL \
      --positions 14,15,16,17,18,19,20,23,48,49,51,52,53,55,56,57,58,59,60,63 \
      --output mutated_sequences.fasta \
      --omit-cys

Notes:
- Positions are 1-based (e.g., 1 is the first residue).
- By default, uses the 20 standard amino acids and skips the wild-type residue
  at each position (yielding 19 mutants per site).
- You can optionally override the amino acid alphabet via --alphabet.
- The --positions argument accepts comma-separated integers and simple ranges
  like "14-20"; both can be mixed (e.g., "14-16,20,23-25").
- Use --omit-cys to exclude cysteine from the mutation alphabet.
"""

import argparse
import sys
from typing import List

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids


def parse_positions(pos_text: str) -> List[int]:
    """Parse a comma-separated positions string that may include ranges."""
    if not pos_text:
        raise ValueError("--positions must be a non-empty string of indices")

    positions = []
    for chunk in pos_text.split(','):
        chunk = chunk.strip()
        if not chunk:
            continue
        if '-' in chunk:
            a, b = chunk.split('-', 1)
            try:
                start = int(a)
                end = int(b)
            except ValueError:
                raise ValueError(f"Invalid range '{chunk}' in --positions")
            if start <= 0 or end <= 0:
                raise ValueError("Positions must be positive (1-based)")
            if start > end:
                start, end = end, start
            positions.extend(range(start, end + 1))
        else:
            try:
                idx = int(chunk)
            except ValueError:
                raise ValueError(f"Invalid index '{chunk}' in --positions")
            if idx <= 0:
                raise ValueError("Positions must be positive (1-based)")
            positions.append(idx)

    seen = set()
    ordered_unique = []
    for p in positions:
        if p not in seen:
            seen.add(p)
            ordered_unique.append(p)
    return ordered_unique


def write_point_mutants(sequence: str, positions: List[int], output_file: str, alphabet: str = STANDARD_AA, omit_cys: bool = False) -> None:
    """Write all single-point mutants (excluding WT residue) to FASTA."""
    if not sequence:
        raise ValueError("Sequence must be a non-empty string")
    if any(ch for ch in sequence if ch.isspace()):
        sequence = "".join(sequence.split())

    seq = sequence.upper()
    aa_set = alphabet.upper()

    if omit_cys:
        aa_set = aa_set.replace("C", "")

    if any(res not in aa_set for res in seq):
        sys.stderr.write(
            "[warning] Sequence contains characters not in the provided alphabet; "
            "those residues can still be mutated using the alphabet provided.\n"
        )

    count = 0
    with open(output_file, 'w') as f:
        for pos in positions:
            index = pos - 1
            if not (0 <= index < len(seq)):
                sys.stderr.write(f"[warning] Position {pos} is out of range for sequence length {len(seq)}. Skipping.\n")
                continue

            original_residue = seq[index]
            for new_residue in aa_set:
                if new_residue == original_residue:
                    continue
                mutated = list(seq)
                mutated[index] = new_residue
                mutated_seq = "".join(mutated)
                header = f">{original_residue}{pos}{new_residue}"
                f.write(f"{header}\n{mutated_seq}\n")
                count += 1

    sys.stderr.write(f"[info] Wrote {count} mutant sequences to '{output_file}'.\n")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate all single point mutants (19/site) for specified positions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sequence", required=True, help="Input sequence (string). Whitespace will be stripped; case-insensitive.")
    p.add_argument("--positions", required=True, help="Comma-separated 1-based positions to mutate. Ranges allowed (e.g., '14-16,20,23-25').")
    p.add_argument("--output", required=True, help="Output FASTA filename.")
    p.add_argument("--alphabet", default=STANDARD_AA, help="Alphabet of residues to mutate to (uppercase string). Defaults to 20 standard amino acids.")
    p.add_argument("--omit-cys", action="store_true", help="Omit cysteine from the mutation alphabet.")
    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        positions = parse_positions(args.positions)
        write_point_mutants(
            sequence=args.sequence,
            positions=positions,
            output_file=args.output,
            alphabet=args.alphabet,
            omit_cys=args.omit_cys,
        )
    except Exception as e:
        sys.stderr.write(f"[error] {e}\n")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
