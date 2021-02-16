#!/usr/bin/env python
import sys

from cupcake.sequence.BioReaders import GMAPSAMReader
from cupcake.sequence.GFF import write_collapseGFF_format


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Convert SAM to collapsed GFF format")
    parser.add_argument("sam_filename")

    args = parser.parse_args()

    if not args.sam_filename.endswith(".sam"):
        print("Only accepts files ending in .sam. Abort!", file=sys.stderr)
        sys.exit(-1)

    prefix = args.sam_filename[:-4]
    output_gff = prefix + ".collapsed.gff"

    with open(output_gff, "w") as f:
        reader = GMAPSAMReader(args.sam_filename, True)
        for r in reader:
            if r.sID == "*":
                continue
            r.strand = r.flag.strand
            r.geneid = r.qID
            r.seqid = r.qID
            r.chr = r.sID
            r.ref_exons = r.segments
            r.start = r.sStart
            r.end = r.sEnd
            r.cds_exons = None
            write_collapseGFF_format(f, r)

    print(f"Output written to {output_gff}.", file=sys.stderr)


if __name__ == "__main__":
    main()
