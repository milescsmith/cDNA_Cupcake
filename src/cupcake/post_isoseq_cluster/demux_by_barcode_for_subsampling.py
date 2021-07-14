#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import re
import sys
from collections import Counter
from csv import DictReader, DictWriter
from typing import Tuple

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger

FIELDNAMES = ["pbid", "pbgene", "length", "refisoform", "refgene", "fl_count"]
pbid_rex = re.compile(r"PB.(\d+).\d+")


app = typer.Typer(name="cupcake.post_isoseq_cluster.demux_by_barcode_for_subsampling")


def demux_for_subsamping(
    class_filename,
    fasta_filename,
    demux_count_file,
    output_prefix,
    out_group_dict,
    ignore_novel,
):
    # read SQANTI classification to get known gene/transcript name
    d = {}  # pbid --> record
    for r in DictReader(open(class_filename), delimiter="\t"):
        d[r["isoform"]] = r

    # get read lengths
    lens = {}  # pbid -> length
    for r in SeqIO.parse(open(fasta_filename), "fasta"):
        lens[r.id] = len(r.seq)

    writers = {}
    handles = {}
    out_groups = set(out_group_dict.values())
    for g in out_groups:
        handles[g] = open(
            f"{output_prefix}_{g}_only.{'ignore_novel' if ignore_novel else 'use_novel'}.for_subsampling.txt",
            "w",
        )
        writers[g] = DictWriter(handles[g], FIELDNAMES, delimiter="\t")
        writers[g].writeheader()

    reader = DictReader(open(demux_count_file), delimiter=",")
    for r in reader:
        if r["id"] not in d:
            logger.info(
                f"WARNING: skipping {r['id']} because not in {class_filename}",
            )
            continue

        m = pbid_rex.match(r["id"])
        if m is None:
            logger.error(
                f"ERROR: unable to parse ID {r['id']}. Expected format PB.X.Y!",
            )
            sys.exit(-1)

        newrec = {"pbid": r["id"], "pbgene": m.group(1), "length": lens[r["id"]]}

        gene = d[r["id"]]["associated_gene"]
        trans = d[r["id"]]["associated_transcript"]
        if gene.startswith("novel") and ignore_novel:
            gene = "NA"
        if trans.startswith("novel"):
            if ignore_novel:
                trans = "NA"
            else:
                trans += r[
                    "id"
                ]  # add an unique identified to make this "novel" refgene unique
        newrec["refgene"] = gene
        newrec["refisoform"] = trans

        group_counts = Counter()
        for b, g in out_group_dict.items():
            group_counts[g] += int(r[b])

        for g in out_groups:
            newrec["fl_count"] = group_counts[g]
            writers[g].writerow(newrec)

    for h in handles.values():
        h.close()


@app.command(name="")
def main(
    class_filename: str = typer.Argument(..., help="SQANTI classification file"),
    fasta_filename: str = typer.Argument(..., help="FASTA filename"),
    demux_count_file: str = typer.Argument(..., help="Demux count file"),
    output_prefix: str = typer.Argument(..., help="Output prefix for GFF outputs"),
    outgroup_dict: Tuple[str, str] = typer.Argument(
        ..., help="Tuples indicating barcode grouping"
    ),
    ignore_novel: bool = typer.Option(
        False,
        help="Ignore novel genes/transcripts (default: off)",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    out_group_dict = dict(eval(outgroup_dict))
    demux_for_subsamping(
        class_filename,
        fasta_filename,
        demux_count_file,
        output_prefix,
        out_group_dict,
        ignore_novel,
    )


if __name__ == "__main__":
    typer.run(main)
