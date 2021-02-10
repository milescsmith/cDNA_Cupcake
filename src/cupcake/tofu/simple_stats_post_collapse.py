#!/usr/bin/env python

import sys
from pathlib import Path
from ..sequence.GFF import collapseGFFReader


def simple_stats_post_collapse(input_prefix):
    input_gff = f"{input_prefix}.gff"

    if not Path(input_gff).exists():
        print(f"Looking for input GFF {input_gff} but not found! Abort!")
        sys.exit(-1)

    with open(f'{input_prefix}.simple_stats.txt', 'w') as f1, open(
        f'{input_prefix}.exon_stats.txt', 'w'
    ) as f2:

        f1.write("pbid\tlocus\tlength\tnum_exon\n")
        f2.write("pbid\texon_index\texon_size\tintron_size\n")

        for r in collapseGFFReader(input_gff):
            f1.write(r.seqid + '\t')
            f1.write(r.seqid.split('.')[1] + '\t')
            sum_len = 0
            for i, e in enumerate(r.ref_exons):
                exon_len = e.end - e.start
                sum_len += exon_len
                f2.write(f"{r.seqid}\t{i+1}\t{exon_len}\t")
                if i == 0:
                    f2.write("NA\n")
                else:
                    f2.write(f"{str(e.start - r.ref_exons[i - 1].end)}\n")

            f1.write(f"{str(sum_len)}\t")
            f1.write(f"{str(len(r.ref_exons))}\n")
        print(f"Output written to: {f1.name},{f2.name}\n")


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("input_prefix", help="Input prefix, ex: hq.5merge.collapsed")

    args = parser.parse_args()
    simple_stats_post_collapse(args.input_prefix)


if __name__ == "__main__":
    main()
