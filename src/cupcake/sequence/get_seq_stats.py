#!/usr/bin/env python

from pathlib import Path

import numpy as np
import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake import cupcake_logger as logger

app = typer.Typer(
    name="cupcake.sequence.get_seq_stats",
    help="Summarize sequence lengths in fasta/fastq",
)


def type_fa_or_fq(filename):
    if filename.suffix.upper() in (".FA", ".FASTA", "CLIPS"):
        return "fasta"
    else:
        return "fastq"


def get_seq_stats(filename, binwidth):
    print("file type is:", type_fa_or_fq(filename))

    with open(f"{filename.name}.seqlengths.txt", "w") as f:
        lens = []
        for r in SeqIO.parse(open(filename), type_fa_or_fq(filename)):
            f.write(f"{r.id}	{str(len(r.seq))}\n")
            lens.append(len(r.seq))

    logger.info(f"{len(lens)} sequences")
    logger.info(f"min: {min(lens)}")
    logger.info(f"max: {max(lens)}")
    logger.info(f"avg: {sum(lens) * 1.0 / len(lens)}")

    # print by 1 kb bins
    logger.info("Length Breakdown by kb range:")

    _max = (max(lens) // binwidth) + 1
    bins = [0] * _max
    for x in lens:
        bins[x // binwidth] += 1

    for i in range(0, _max):
        if binwidth == 1000:
            print(f"{i}-{i + 1} kb: {bins[i]}")
        else:
            print(f"{i * binwidth}-{(i + 1) * binwidth}: {bins[i]}")

    print("5-95% percentile:", np.percentile(lens, 5), np.percentile(lens, 95))


@app.command(name="")
def main(
    filename: str = typer.Argument(..., help="Input fasta/fastq filename"),
    binwidth: int = typer.Option(
        1000,
        "--binwidth",
        "-b",
        help="Bin width, in bp (default: 1000 bp)",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:

    get_seq_stats(Path(filename), binwidth)


if __name__ == "__main__":
    typer.run(main)
