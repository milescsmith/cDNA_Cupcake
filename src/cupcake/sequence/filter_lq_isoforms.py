#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Sequence headers must have:

@i0_LQ_sample32f89e|c31/f12p0/433 isoform=c31;full_length_coverage=12;non_full_length_coverage=0;isoform_length=433;expected_accuracy=0.300

"""
import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger

app = typer.Typer(name="cupcake.sequence.filter_lq_isoforms")


def func(r):
    fl_count = None
    exp_acc = None
    for x in r.description.split(";"):
        if x.find("=") == -1:
            continue
        a, b = x.split("=")
        if a == "full_length_coverage":
            fl_count = int(b)
        elif a == "expected_accuracy":
            exp_acc = float(b)
    return fl_count, exp_acc


def filter_lq_isoforms(
    fastq_filename: str,
    output_filename: str,
    min_fl_count: int,
    min_exp_acc: float,
    is_flnc: bool,
) -> None:
    for r in SeqIO.parse(open(fastq_filename), "fastq"):
        fl_count, exp_acc = func(r)
        if not is_flnc and fl_count is None:
            raise RuntimeError(
                "Sequence header does not include field `full_length_coverage=`. Abort!"
            )
        if exp_acc is None:
            raise RuntimeError(
                "Sequence header does not include field `expected_accuracy=`. Please run calc_expected_accuracy_from_fastq.py script first!"
            )
        if (is_flnc or fl_count >= min_fl_count) and exp_acc >= min_exp_acc:
            logger.info(
                f"Including {r.id} into output file (FL: {fl_count}, acc: {exp_acc})."
            )
            SeqIO.write(r, output_filename, "fastq")


@app.command(name="")
def main(
    fastq_filename: str = typer.Argument(
        ..., help="LQ FASTQ filename (ex: lq_isoforms.fastq)"
    ),
    output_filename: str = typer.Argument(..., help="Output FASTQ filename"),
    min_fl_count: int = typer.Option(2, help="Minimum FL count (default: 2)."),
    min_exp_acc: float = typer.Option(
        0.99, help="Minimum predicted accuracy (default: 0.99)."
    ),
    is_flnc: bool = typer.Option(False, help="Input FASTQ is FLNC, not LQ"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    if not 0 <= min_exp_acc <= 1:
        raise ValueError("min_exp_acc much be between 0 and 1")

    filter_lq_isoforms(
        fastq_filename,
        output_filename,
        min_fl_count,
        min_exp_acc,
        is_flnc,
    )


if __name__ == "__main__":
    typer.run(main)
