#!/usr/bin/env python

import random

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake import cupcake_logger as logger

app = typer.Typer(
    name="cupcake.sequence.randomly_select_sequences",
    help="Randomly select N sequences from fasta/fastq files",
)


def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith(".FA") or file.endswith(".FASTA"):
        return "fasta"
    else:
        return "fastq"


def sep_by_primer(filename, output_prefix, sample_size):
    filetype = type_fa_or_fq(filename)

    ids = [r.id for r in SeqIO.parse(open(filename), filetype)]

    n = len(ids)
    if sample_size > n:
        logger.warning(
            f"WARNING: {filename} contains only {n} sequences but subsample size at {sample_size}! Simply output whole file."
        )

    chosen_ids = random.sample(ids, min(n, sample_size))

    with open(f"{output_prefix}.random{str(sample_size)}.{filetype}", "w") as f:
        for r in SeqIO.parse(open(filename), filetype):
            if r.id in chosen_ids:
                SeqIO.write(r, f, filetype)

        logger.info(f"Randomly selected sequences written to {f.name}.")


@app.command(name="")
def main(
    filename: str = typer.Argument(..., help="Input fasta/fastq filename"),
    output_prefix: str = typer.Argument(..., help="Output file prefix"),
    sample_size: int = typer.Argument(..., help="Subsample size"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:

    sep_by_primer(filename, output_prefix, sample_size)


if __name__ == "__main__":
    typer.run(main)
