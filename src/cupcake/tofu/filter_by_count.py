#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import os
import sys
from csv import DictReader, DictWriter

import typer
from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import GFF

"""
Given the collapse script result, further filter by FL counts.

Input: <prefix> (must have .group.txt, .abundance.txt, .gff, .rep.fq)
Output: <output>.min_fl_<threshold>.{group|abundance|gff|rep fq}
"""


app = typer.Typer(
    name="filter_by_count",
    help="Given the collapse script result, further filter by FL counts.",
)


def filter_by_count(
    input_prefix: str,
    output_prefix: str,
    min_count: int,
    dun_use_group_count: bool = False,
):

    group_filename = f"{input_prefix}.group.txt"
    count_filename = f"{input_prefix}.abundance.txt"
    gff_filename = f"{input_prefix}.gff"
    rep_filenames = [
        (f"{input_prefix}.rep.fq", "fastq"),
        (f"{input_prefix}.rep.fastq", "fastq"),
        (f"{input_prefix}.rep.fa", "fasta"),
        (f"{input_prefix}.rep.fasta", "fasta"),
    ]

    rep_filename = None
    rep_type = None
    for x, feature in rep_filenames:
        if os.path.exists(x):
            rep_filename = x
            rep_type = feature

    if rep_filename is None:
        logger.error(
            f"Expected to find input fasta or fastq files {input_prefix}.rep.fa or {input_prefix}.rep.fq. Not found. Abort!"
        )
        sys.exit(-1)

    if not dun_use_group_count:
        # read group
        group_max_count_fl = {}
        group_max_count_p = {}
        for line in open(group_filename):
            # ex: PB.1.1  i0HQ_54b0ca|c58773/f30p16/700
            pbid, members = line.strip().split("\t")
            group_max_count_fl[pbid] = 0
            group_max_count_p[pbid] = 0
            members = members.split(",")
            for m in members:
                i = m.find("|")
                if i > 0:
                    tmp = m.split("|")[1].split("/")[1]  # ex: tmp = f30p16
                else:
                    tmp = m.split("/")[1]
                fl_count, p_count = tmp.split("p")
                fl_count = int(fl_count[1:])
                p_count = int(p_count)
                group_max_count_fl[pbid] = max(group_max_count_fl[pbid], fl_count)
                group_max_count_p[pbid] = max(group_max_count_p[pbid], p_count)

    # read abundance first
    with open(count_filename) as f:
        count_header = ""
        while True:
            cur_pos = f.tell()
            line = f.readline()
            if not line.startswith("#"):
                f.seek(cur_pos)
                break
            else:
                count_header += line
        d = {r["pbid"]: r for r in DictReader(f, delimiter="\t")}
        for k, v in d.items():
            print(k, v)

    # group_max_count_p NOT used for now
    good = [
        x
        for x in d
        if int(d[x]["count_fl"]) >= min_count
        and (dun_use_group_count or group_max_count_fl[x] >= min_count)
    ]

    # write output GFF
    with open(f"{output_prefix}.gff", "w") as f:
        for r in GFF.collapseGFFReader(gff_filename):
            if r.seqid in good:
                GFF.write_collapseGFF_format(f, r)

    # write output rep.fq
    with open(
        f"{output_prefix}.rep.{('fq' if rep_type == 'fastq' else 'fa')}", "w"
    ) as f:
        for r in SeqIO.parse(open(rep_filename), rep_type):
            if r.name.split("|")[0] in good:
                SeqIO.write(r, f, rep_type)

    # write output to .abundance.txt
    with open("{output_prefix}.abundance.txt", "w") as f:
        f.write(count_header)
        writer = DictWriter(
            f,
            fieldnames=[
                "pbid",
                "count_fl",
                "count_nfl",
                "count_nfl_amb",
                "norm_fl",
                "norm_nfl",
                "norm_nfl_amb",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for k in good:
            r = d[k]
            writer.writerow(r)

    logger.info(
        f"Output written to: {output_prefix}.gff\n"
        f"Output written to: {rep_filename}\n"
        f"Output written to: {output_prefix}.abundance.txt"
    )


@app.command(name="")
def main(
    input_prefix: str = typer.Argument(...),
    min_count: int = typer.Optional(2, help="Minimum FL count (default: 2)"),
    dun_use_group_count: bool = typer.Optional(
        True,
        "--dun_use-group_count",
        help="Turn off more stringent min count (default: off)",
    ),
) -> None:

    output_prefix = f"{input_prefix}.min_fl_{min_count}"
    filter_by_count(input_prefix, output_prefix, min_count, dun_use_group_count)


if __name__ == "__main__":
    typer.run(main)
