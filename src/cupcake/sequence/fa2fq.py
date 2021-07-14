#!/usr/bin/env python
from pathlib import Path

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake import cupcake_logger as logger

app = typer.Typer(name="cupcake.sequence.fa2fq", help="Convert fasta to fastq")


def fa2fq(input_file):
    if not (
        input_file.lower().endswith(".fasta") or input_file.lower().endswith(".fa")
    ):
        raise AssertionError(
            f"Input {input_file} does not end with .fasta or .fa! Abort"
        )
    output = Path(input_file).with_suffix(".fastq")

    f = open(output, "w")
    for r in SeqIO.parse(open(input_file), "fasta"):
        r.letter_annotations["phred_quality"] = [60] * len(r.seq)
        SeqIO.write(r, f, "fastq")
    f.close()

    logger.info(f"Output written to {f.name}")
    return f.name


@app.command(name="")
def main(
    fasta_filename: str = typer.Argument(
        ..., help="input fasta (must end with .fasta or .fa)"
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    fa2fq(fasta_filename)


if __name__ == "__main__":
    typer.run(main)
