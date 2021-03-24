#!/usr/bin/env python

__version__ = "1.0"

import sys

import typer
from Bio import SeqIO
from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.sequence.fq2fa", help="Convert fastq to fasta")


def fq2fa(input):
    if not input.lower().endswith(".fastq") or input.lower().endswith(".fq"):
        raise AssertionError(f"Input {input} does not end with .fastq or .fq! Abort")

    output = f"{input[:input.rfind('.')]}.fasta"

    for r in SeqIO.parse(open(input), "fastq"):
        SeqIO.write(r, output, "fasta")

    logger.info(f"Output written to {output}")


@app.command(name="")
def main(
    fastq_filename: str = typer.Argument(
        ..., help="input fastq (must end with .fastq or .fq)"
    )
):
    fq2fa(fastq_filename)


if __name__ == "__main__":
    typer.run(main)
