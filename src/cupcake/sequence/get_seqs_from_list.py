#!/usr/bin/env python

import sys
from pathlib import Path

import typer
from Bio import SeqIO

app = typer.Typer(
    name="cupcake.sequence.get_seqs_from_list",
    help="Get sequences from a fasta/fastq file from a list",
)


def get_seqs_from_list(fafq, listfile, partial_ok=False, exclude=False):
    f = sys.stdout
    fafq = Path(fafq)
    filetype = "fastq" if fafq.suffix.upper() in (".FQ", ".FASTQ") else "fasta"
    seqs = {line.strip() for line in open(listfile)}
    for r in SeqIO.parse(open(fafq), filetype):
        id_seen = (
            r.id in seqs
            or r.id.split("|")[0] in seqs
            or (partial_ok and any(r.id.startswith(x) for x in seqs))
        )
        if id_seen ^ exclude:
            SeqIO.write(r, f, filetype)


@app.command(name="")
def main(
    fasta_filename: str = typer.Argument(
        ..., help="Input fasta/fastq filename to extract sequences from"
    ),
    list_filename: str = typer.Argument(..., help="List of sequence IDs to extract"),
    partial: bool = typer.Option(
        False,
        help="OK if seq IDs only match the beginning",
    ),
    exclude: bool = typer.Option(
        False,
        help="Output sequences NOT in the list, default OFF",
    ),
) -> None:

    get_seqs_from_list(fasta_filename, list_filename, partial, exclude)


if __name__ == "__main__":
    typer.run(main)
