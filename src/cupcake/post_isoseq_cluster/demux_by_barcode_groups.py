#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Given a pooled input GFF + demux CSV file, write out per-{barcode group} GFFs
If input fasta/fastq is given, optionally also output per-{barcode group} FASTA/FASTQ
"""
import re
from collections import defaultdict
from csv import DictReader
from typing import Optional

import cupcake.sequence.GFF as GFF
import typer
from Bio import SeqIO

rex_pbid = re.compile(r"(PB.\d+.\d+)(|\S+)")


app = typer.Typer(name="cupcake.post_isoseq_cluster.demux_by_barcode_groups")


def get_type_fafq(in_filename):
    in_filename = in_filename.upper()
    if in_filename.endswith(".FA") or in_filename.endswith("FASTA"):
        return "fasta"
    elif in_filename.endswith(".FQ") or in_filename.endswith("FASTQ"):
        return "fastq"
    else:
        raise Exception(
            f"Unrecognized file suffix .{in_filename[in_filename.find('.'):]}! Must end with .fasta or .fastq!"
        )


def regroup_gff(
    pooled_gff, demux_count_file, output_prefix, out_group_dict, in_fafq=None
):
    """
    :param pooled_sam: SAM file
    :param demux_count_file: comma-delimited per-barcode count file
    :param output_prefix: output prefix for GFF
    :param out_group_dict: dict of barcode name --> group to be long in  (ex: {'EM1':'EM', 'EM2':'EM'})
    :param in_fafq: optional fasta/fastq that was input to SAM
    """
    if in_fafq is not None:
        type_fafq = get_type_fafq(in_fafq)
    in_tissue = defaultdict(
        lambda: set()
    )  # pbid --> list of tissue it is in (EM, END, R)

    for r in DictReader(open(demux_count_file), delimiter=","):
        for k, v in r.items():
            if k != "id" and int(v) > 0:
                in_tissue[r["id"]].add(k)

    # in_tissue = dict(in_tissue)

    handles = {}
    handles_fafq = {}
    for g in out_group_dict.values():
        handles[g] = open(f"{output_prefix}_{g}_only.gff", "w")
        if in_fafq is not None:
            handles_fafq[g] = open(f"{output_prefix}_{g}_only.{type_fafq}", "w")

    if in_fafq is not None:
        fafq_dict = SeqIO.to_dict(SeqIO.parse(open(in_fafq), type_fafq))
        fafq_dict_keys = list(fafq_dict.keys())
        for k in fafq_dict_keys:
            m = rex_pbid.match(k)
            if m is not None:
                fafq_dict[m.group(1)] = fafq_dict[k]
    reader = GFF.collapseGFFReader(pooled_gff)
    for r in reader:
        groups_to_write_in = set()
        pbid = r.seqid
        if pbid not in in_tissue:
            logger.info(
                f"WARNING: {pbid} does not belong to any group indicated by outgroup_dict"
            )
        for tissue in in_tissue[pbid]:
            groups_to_write_in.add(out_group_dict[tissue])

        for g in groups_to_write_in:
            GFF.write_collapseGFF_format(handles[g], r)
            if in_fafq is not None:
                SeqIO.write(fafq_dict[pbid], handles_fafq[g], type_fafq)


@app.command(name="")
def main(
    pooled_gff: str = typer.Argument(..., help="Pooled GFF file"),
    demux_count_file: str = typer.Argument(..., help="Demux count file"),
    output_prefix: str = typer.Argument(..., help="Output prefix for GFF outputs"),
    outgroup_dict: str = typer.Argument(..., help="Tuples indicating barcode grouping"),
    pooled_fastx: Optional[str] = typer.Option(
        None,
        help="Pooled FASTA/FASTQ (optional, if given, will also output demux fa/fq)",
    ),
) -> None:
    tmp = eval(outgroup_dict)
    out_group_dict = dict([tmp]) if len(tmp) == 1 else dict(tmp)
    regroup_gff(
        pooled_gff,
        demux_count_file,
        output_prefix,
        out_group_dict,
        pooled_fastx,
    )


if __name__ == "__main__":
    main()
