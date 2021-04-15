#!/usr/bin/env python
from pathlib import Path

import typer
from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.sequence.fq2fa", help="Convert fastq to fasta")


def fq2fa(input_file):
    if not input_file.lower().endswith(".fastq") or input_file.lower().endswith(".fq"):
        raise AssertionError(
            f"Input {input_file} does not end with .fastq or .fq! Abort"
        )

    output = Path(input_file).with_suffix(".fasta")

    for r in SeqIO.parse(open(input_file), "fastq"):
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
