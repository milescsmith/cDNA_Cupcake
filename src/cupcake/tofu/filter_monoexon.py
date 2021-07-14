#!/usr/bin/env python

__author__ = "etseng@pacb.com"

"""
Filter away transcripts that are mono-exonic.

Required input: <input_prefix>.gff, .rep.fq, .abundance.txt

Example:
    filter_monoexon.py test.collapsed test.collapsed.nomono

"""

import sys
from csv import DictReader, DictWriter
from pathlib import Path
from typing import Dict, Tuple

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger
from cupcake.sequence import GFF

app = typer.Typer(name="cupcake.tofu.filter_monoexon")


def sanity_check_collapse_input(input_prefix: str) -> Tuple[Path, Path, Path]:
    """
    Check that
    1. the count, gff, rep files exist
    2. the number of records agree among the three
    """
    # group_filename =  f"{input_prefix}.group.txt"
    count_filename = Path(f"{input_prefix}.abundance.txt")
    gff_filename = Path(f"{input_prefix}.gff")
    rep_filename = Path(f"{input_prefix}.rep.fq")
    if not count_filename.exists():
        logger.error(f"File {count_filename} does not exist. Abort!")
        sys.exit(-1)
    if not gff_filename.exists():
        logger.error(f"File {gff_filename} does not exist. Abort!")
        sys.exit(-1)
    if not rep_filename.exists():
        logger.error(f"File {rep_filename} does not exist. Abort!")
        sys.exit(-1)

    pbids1 = {[r.id for r in SeqIO.parse(open(rep_filename, "r"), "fastq")]}
    pbids2 = {[r.seqid for r in GFF.collapseGFFReader(gff_filename)]}
    pbids3 = {read_count_file(count_filename)[0].keys()}

    if (
        len(pbids1) != len(pbids2)
        or len(pbids2) != len(pbids3)
        or len(pbids1) != len(pbids3)
    ):
        logger.error(
            "The number of PBID records in the files disagree! Sanity check failed."
        )
        logger.error(f"# of PBIDs in {rep_filename}: {len(pbids1)}")
        logger.error(f"# of PBIDs in {gff_filename}: {len(pbids2)}")
        logger.error(f"# of PBIDs in {count_filename}: {len(pbids3)}")
        sys.exit(-1)

    return count_filename, gff_filename, rep_filename


def read_count_file(count_filename: Path) -> Tuple[Dict[str, str], str]:
    with open(count_filename, "r") as f:
        count_header = ""
        while True:
            cur_pos = f.tell()
            line = f.readline()
            if not line.startswith("#"):
                f.seek(cur_pos)
                break
            else:
                count_header += line
        d = {r["pbid"]: r for r in DictReader(f, delimiter="\t")}

    return d, count_header


@app.command(name="")
def main(
    input_prefix: str = typer.Argument(
        ..., help="Input prefix (ex: test.collapsed.min_fl_2)"
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:

    output_prefix = f"{input_prefix}.nomono"

    count_filename, gff_filename, rep_filename = sanity_check_collapse_input(
        input_prefix
    )

    good = []
    with open(f"{output_prefix}.gff", "w") as f:
        reader = GFF.collapseGFFReader(gff_filename)
        for r in reader:
            assert r.seqid.startswith("PB.")
            if len(r.ref_exons) > 1:
                good.append(r.seqid)
                GFF.write_collapseGFF_format(f, r)

    # read abundance first
    d, count_header = read_count_file(count_filename)

    # write output rep.fq
    with open(f"{output_prefix}.rep.fq", "w") as f:
        for r in SeqIO.parse(open(rep_filename, "r"), "fastq"):
            if r.name.split("|")[0] in good:
                SeqIO.write(r, f, "fastq")

    # write output to .abundance.txt
    with open(f"{output_prefix}.abundance.txt", "w") as f:
        f.write(count_header)
        writer = DictWriter(
            f,
            fieldnames=[
                "pbid",
                "count_fl",
                "count_nfl",
                "count_nfl_amb",
                "norm_fl",
                "norm_nfl",
                "norm_nfl_amb",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for k in good:
            r = d[k]
            writer.writerow(r)

    logger.info(f"Output written to:{output_prefix}.gff")
    logger.info(f"Output written to:{output_prefix}.rep.fq")


if __name__ == "__main__":
    typer.run(main)
