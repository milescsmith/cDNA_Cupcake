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

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake import cupcake_logger as logger

rex_pbid = re.compile(r"(PB.\d+).(\d+)")


app = typer.Typer(
    name="cupcake.sequence.group_ORF_sequences", help="De-duplicate ORF FAA file."
)


def dedup_ORFs(faa_filename, output_prefix, is_pbid):

    seq_dict = OrderedDict()  # ORF seq --> list of IDs

    pbid_counter = Counter()  # PB.X  --> counter  (not used if is_pbid is False)

    for r in SeqIO.parse(open(faa_filename), "fasta"):
        s = str(r.seq).upper()
        if s not in seq_dict:
            seq_dict[s] = []
        seq_dict[s].append(r.id)

    with open(f"{output_prefix}.faa", "w") as f1, open(
        f"{output_prefix}.group.txt", "w"
    ) as f2:

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


@app.command(name="")
def main(
    input_faa: str = typer.Argument(..., help="Input FAA filename"),
    output_prefix: str = typer.Argument(..., help="Output prefix"),
    is_pbid: bool = typer.Option(
        False,
        help="FAA IDs are in PB.X.Y format (default: off)",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:

    dedup_ORFs(input_faa, output_prefix, is_pbid)


if __name__ == "__main__":
    typer.run(main)
