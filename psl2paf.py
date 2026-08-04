#!/usr/bin/env python3
"""
Convert halSynteny PSL output to standard PAF format.

PSL columns (0-indexed):
  0  matches
  1  misMatches
  2  repMatches
  3  nCount
  4  qNumInsert
  5  qBaseInsert
  6  tNumInsert
  7  tBaseInsert
  8  strand
  9  qName
  10 qSize
  11 qStart
  12 qEnd
  13 tName
  14 tSize
  15 tStart
  16 tEnd
  17 blockCount
  18 blockSizes
  19 qStarts
  20 tStarts

PAF columns (0-indexed):
  0  qName
  1  qLen
  2  qStart
  3  qEnd
  4  strand
  5  tName
  6  tLen
  7  tStart
  8  tEnd
  9  numMatches
  10 alignLen
  11 mapq (set to 255 = missing)
"""

import sys
import argparse

def psl_to_paf(infile, outfile, min_match=0):
    for line in infile:
        line = line.rstrip("\n")
        if not line or line.startswith("psLayout") or line.startswith("match") or line.startswith("---") or line.startswith("#"):
            continue

        f = line.split("\t")
        if len(f) < 21:
            continue

        try:
            matches     = int(f[0])
            mismatches  = int(f[1])
            rep_matches = int(f[2])
            strand_raw  = f[8]
            q_name      = f[9]
            q_size      = int(f[10])
            q_start     = int(f[11])
            q_end       = int(f[12])
            t_name      = f[13]
            t_size      = int(f[14])
            t_start     = int(f[15])
            t_end       = int(f[16])
            block_count = int(f[17])
        except (ValueError, IndexError):
            continue

        if matches < min_match:
            continue

        # PSL strand field can be "+", "-", "++", "+-", "-+", "--"
        # First char = query strand, second (if present) = target strand
        # PAF strand = relative orientation between query and target
        if len(strand_raw) == 2:
            q_strand = strand_raw[0]
            t_strand = strand_raw[1]
            # If target strand is '-', coordinates are on reverse complement
            # and we need to flip query strand for PAF convention
            if t_strand == '-':
                strand = '-' if q_strand == '+' else '+'
                # Reverse-complement target coords
                t_start, t_end = t_size - t_end, t_size - t_start
            else:
                strand = q_strand
        else:
            strand = strand_raw if strand_raw in ('+', '-') else '+'

        align_len = matches + mismatches + rep_matches
        num_matches = matches + rep_matches

        # mapq: not available from PSL, use 255 (missing)
        mapq = 255

        out = [
            q_name,
            str(q_size),
            str(q_start),
            str(q_end),
            strand,
            t_name,
            str(t_size),
            str(t_start),
            str(t_end),
            str(num_matches),
            str(align_len),
            str(mapq),
        ]
        outfile.write("\t".join(out) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert halSynteny PSL output to PAF format."
    )
    parser.add_argument(
        "input", nargs="?", default="-",
        help="Input PSL file (default: stdin)"
    )
    parser.add_argument(
        "output", nargs="?", default="-",
        help="Output PAF file (default: stdout)"
    )
    parser.add_argument(
        "--min-match", type=int, default=0,
        help="Minimum number of matching bases to include (default: 0)"
    )
    args = parser.parse_args()

    infile  = open(args.input,  "r") if args.input  != "-" else sys.stdin
    outfile = open(args.output, "w") if args.output != "-" else sys.stdout

    try:
        psl_to_paf(infile, outfile, min_match=args.min_match)
    finally:
        if args.input  != "-": infile.close()
        if args.output != "-": outfile.close()


if __name__ == "__main__":
    main()
