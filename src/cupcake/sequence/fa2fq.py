#!/usr/bin/env python

__version__ = "1.0"

import typer
from Bio import SeqIO
from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.sequence.fa2fq", help="Convert fasta to fastq")


def fa2fq(input):
    if not input.lower().endswith(".fasta") or input.lower().endswith(".fa"):
        raise AssertionError(f"Input {input} does not end with .fasta or .fa! Abort")
    output = f"{input[:input.rfind('.')]}.fastq"

    f = open(output, "w")
    for r in SeqIO.parse(open(input), "fasta"):
        r.letter_annotations["phred_quality"] = [60] * len(r.seq)
        SeqIO.write(r, f, "fastq")
    f.close()

    logger.info(f"Output written to {f.name}")
    return f.name


@app.command(name="")
def main(
    fasta_filename: str = typer.Argument(
        ..., help="input fasta (must end with .fasta or .fa)"
    )
) -> None:
    fa2fq(fasta_filename)


if __name__ == "__main__":
    typer.run(main)
