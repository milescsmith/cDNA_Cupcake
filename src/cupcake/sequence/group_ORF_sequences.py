#!/usr/bin/env python
"""
Given a FAA file, group identical ORFs.

INPUT: a single FAA file
OUTPUT: a de-duped FAA file with a companion "group" file

The de-duped FAA format:

>ORFgroup_<PB.X>_<index>  if it is in PB.X.Y format

>ORFgroup_<index>  if not PB.X.Y format

The group file format:

ORFgroup_<index> \t comma-sep list of IDs with the same ORF
"""

import re
from collections import Counter, OrderedDict

from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger

rex_pbid = re.compile(r"(PB.\d+).(\d+)")


def dedup_ORFs(faa_filename, output_prefix, is_pbid):

    seq_dict = OrderedDict()  # ORF seq --> list of IDs

    pbid_counter = Counter()  # PB.X  --> counter  (not used if is_pbid is False)

    for r in SeqIO.parse(open(faa_filename), "fasta"):
        s = str(r.seq).upper()
        if s not in seq_dict:
            seq_dict[s] = []
        seq_dict[s].append(r.id)

    with open(f"{output_prefix}.faa", "w") as f1, open(f"{output_prefix}.group.txt", "w") as f2:

        for i, s in enumerate(seq_dict):
            newid = None
            if is_pbid:
                m = rex_pbid.match(
                    seq_dict[s][0]
                )  # we will just take the first member and use the PB.X.Y
                if m is None:
                    logger.warning(
                        f"WARNING: seqid {seq_dict[s][0]} is not in PB.X.Y format!",
                    )
                else:
                    pb_x = m.group(1)  # ex: PB.10
                    pbid_counter[pb_x] += 1
                    newid = f"ORFgroup_{pb_x}_{pbid_counter[pb_x]}"
            if newid is None:
                newid = f"ORFgroup_{i + 1}"
            f1.write(f">{newid}\n{s}\n")
            f2.write(f"{newid}\t{','.join(seq_dict[s])}\n")

    logger.info(f"Output written to: {f1.name},{f2.name}")


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("De-duplicate ORF FAA file.")
    parser.add_argument("input_faa", help="Input FAA filename")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument(
        "--is_pbid",
        action="store_true",
        default=False,
        help="FAA IDs are in PB.X.Y format (default: off)",
    )

    args = parser.parse_args()

    dedup_ORFs(args.input_faa, args.output_prefix, args.is_pbid)


if __name__ == "__main__":
    main()
