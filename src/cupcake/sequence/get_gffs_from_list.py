#!/usr/bin/env python
import sys

import typer

from cupcake.sequence import GFF

app = typer.Typer(
    name="cupcake.sequence.get_gffs_from_list",
    help="Get records from a GFF file from a list",
)


def get_gff_from_list(gff_filename, listfile, partial_ok=False):
    seqs = [line.strip() for line in open(listfile)]
    for r in GFF.collapseGFFReader(gff_filename):
        if (
            r.seqid in seqs
            or r.seqid.split("|")[0] in seqs
            or (partial_ok and any(r.seqid.startswith(x) for x in seqs))
        ):
            GFF.write_collapseGFF_format(sys.stdout, r)


@app.command(name="")
def main(
    gff_filename: str = typer.Argument(
        ..., help="Input gff filename to extract sequences from"
    ),
    list_filename: str = typer.Argument(..., help="List of sequence IDs to extract"),
    partial: bool = typer.Option(
        False,
        help="OK if seq IDs only match the beginning",
    ),
) -> None:

    get_gff_from_list(gff_filename, list_filename, partial)


if __name__ == "__main__":
    typer.run(main)
