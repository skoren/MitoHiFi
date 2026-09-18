#!/usr/bin/env python3
"""Extract unique hifi_m read IDs for a scaffold using scfmap and layout files."""

import argparse
from pathlib import Path
import re
import sys


def extract(scaffold, scfmap, layout):
    pieces = set()
    selected = False
    found = False
    with open(scfmap) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "path":
                selected = len(fields) > 1 and fields[1] == scaffold
                found |= selected
            elif fields[0] == "end":
                selected = False
            elif selected:
                # Scaffold gaps such as [N100000N] have no layout records.
                pieces.update(field for field in fields
                              if re.fullmatch(r"piece\d+", field))
    if not found:
        raise ValueError("Scaffold {!r} not found in scfmap".format(scaffold))
    if not pieces:
        raise ValueError("Scaffold {!r} has no pieces".format(scaffold))

    seen_pieces = set()
    seen_reads = set()
    reads = []
    selected = False
    with open(layout) as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "tig":
                if len(fields) < 2:
                    raise ValueError("Layout contains a tig line without a piece ID")
                selected = fields[1] in pieces
                if selected:
                    seen_pieces.add(fields[1])
            elif fields[0] == "end":
                selected = False
            elif selected and fields[0].startswith("hifi_m"):
                read = fields[0]
                if read not in seen_reads:
                    reads.append(read.removeprefix("hifi_"))
                    seen_reads.add(read)
    missing = pieces - seen_pieces
    if missing:
        raise ValueError("Pieces missing from layout: " + ", ".join(sorted(missing)))
    return reads, len(pieces)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scaffold", help="Scaffold ID, e.g. unassigned-0000497")
    parser.add_argument("scfmap", type=Path)
    parser.add_argument("layout", type=Path)
    parser.add_argument("output", type=Path, help="Output text file (overwritten)")
    args = parser.parse_args()
    if args.output.resolve() in {args.scfmap.resolve(), args.layout.resolve()}:
        parser.error("Output must not overwrite an input file")
    try:
        reads, count = extract(args.scaffold, args.scfmap, args.layout)
        with open(args.output, "a") as handle:
            for read in reads:
                handle.write(read + "\n")
    except (OSError, ValueError) as error:
        parser.exit(1, "Error: {}\n".format(error))
    print("Wrote {} unique HiFi read IDs from {} pieces to {}".format(
        len(reads), count, args.output), file=sys.stderr)


if __name__ == "__main__":
    main()
