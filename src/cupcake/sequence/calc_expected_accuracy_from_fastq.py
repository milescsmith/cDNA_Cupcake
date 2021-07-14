#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Helper script for calculating expected accuracy from FASTQ sequences.

ex:
i0_LQ_sample32f89e|c32/f1p0/929 isoform=c32;full_length_coverage=1;non_full_length_coverage=0;isoform_length=929;expected_accuracy=1.0
"""

import typer
from Bio import SeqIO

from cupcake import version_callback

app = typer.Typer(name="cupcake.sequence.calc_expected_accuracy_from_fastq")


def phred_to_qv(phred):
    """Phred value to quality value."""
    return 10 ** -(phred / 10.0)


def calc_exp_acc(r, qv_trim_5, qv_trim_3):
    """
    :param r: SeqIO Fastq record
    :param qv_trim_5: 5' trimming
    :param qv_trim_3: 3' trimming
    """
    assert qv_trim_5 >= 0 and qv_trim_3 >= 0
    qv = r.letter_annotations["phred_quality"]
    qv_len = len(qv)
    q = [phred_to_qv(x) for x in qv]
    if qv_trim_3 == 0:
        err_sum = sum(q[qv_trim_5:])
    else:
        err_sum = sum(q[qv_trim_5:-qv_trim_3])
    return 1.0 - (err_sum / float(qv_len))


def cal_expected_accuracy_from_fastq(
    fastq_filename, output_filename, qv_trim_5, qv_trim_3
):
    with open(output_filename, "w") as f:
        for r in SeqIO.parse(open(fastq_filename), "fastq"):
            exp_acc = calc_exp_acc(r, qv_trim_5, qv_trim_3)
            r.description += f";expected_accuracy={exp_acc:.3f}"
            SeqIO.write(r, f, "fastq")


@app.command(name="")
def main(
    fastq_filename: str = typer.Argument(
        ..., help="FASTQ filename (ex: lq_isoforms.fastq)"
    ),
    output_filename: str = typer.Argument(..., help="Output FASTQ filename"),
    qv_trim_5: int = typer.Option(
        100,
        help="Ignore length on 5' for QV calculation (default: 100 bp)",
    ),
    qv_trim_3: int = typer.Option(
        30,
        help="Ignore length on 3' for QV calculation (default: 30 bp)",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    cal_expected_accuracy_from_fastq(
        fastq_filename, output_filename, qv_trim_5, qv_trim_3
    )


if __name__ == "__main__":
    typer.run(main)
