#!/usr/bin/env python
from collections import defaultdict

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger

app = typer.Typer(
    name="cupcake.annotation.parse_matchAnnot", help="Parse MatchAnnot result"
)


def type_fa_or_fq(file):
    file = file.upper()
    if file.endswith(".FA") or file.endswith(".FASTA"):
        return "fasta"
    else:
        return "fastq"


def parse_matchAnnot(fa_or_fq, filename, not_pbid=False, parse_FL_coverage=False):
    pbids = []
    fl_cov = {}  # only used if parse_FL_coverage is True
    for r in SeqIO.parse(open(fa_or_fq), type_fa_or_fq(fa_or_fq)):
        _id = r.id if not_pbid else r.id.split("|")[0]
        pbids.append(_id)
        if parse_FL_coverage:
            try:
                cov = int(r.description.split("full_length_coverage=")[1].split(";")[0])
                fl_cov[_id] = cov
            except:
                logger.error(
                    f"WARNING: Unable to extract `full_length_coverage=` from {r.description}. Mark as NA."
                )
                fl_cov[_id] = "NA"

    match = defaultdict(lambda: (None, None, 0))  # ex: PB.1.1 -> (NOC2L, NOC2L-001, 5)

    for line in open(filename):
        i = line.find("result:")
        if i >= 0:
            raw = line[i:].strip().split()
            if len(raw) < 7:
                continue
            pbid = raw[1] if not_pbid else raw[1].split("|")[0]
            gene = raw[2]
            isoform = raw[3]
            score = int(raw[7])
            if score > match[pbid][1]:
                match[pbid] = (gene, isoform, score)

    f = open(f"{filename}.parsed.txt", "w")
    f.write("pbid\tpbgene\trefisoform\trefgene\tscore")
    if parse_FL_coverage:
        f.write("\tcount_fl")
    f.write("\n")
    for pbid in pbids:
        if not_pbid:
            pbpre = pbid
        else:
            pbpre = pbid.split(".")[1]
        _cov_text = f"\t{fl_cov[pbid]}" if parse_FL_coverage else ""
        if pbid not in match:
            f.write(f"{pbid}\t{pbpre}\tNA\tNA\tNA{_cov_text}\n")
        else:
            gene, isoform, score = match[pbid]
            if gene is None:
                f.write(f"{pbid}\t{pbpre}\tNA\tNA\tNA{_cov_text}\n")
            else:
                f.write(f"{pbid}\t{pbpre}\t{isoform}\t{gene}\t{score}{_cov_text}\n")
    f.close()
    logger.info(f"Output written to: {f.name}")


@app.command(name="")
def main(
    fa_or_fq: str = typer.Argument(
        ...,
        help="Fasta/Fastq filename used to create the SAM file for matchAnnot",
    ),
    match_filename: str = typer.Argument(..., help="MatchAnnot filename"),
    not_pbid: bool = typer.Option(
        False,
        help="Turn this on if not sequence ID is not PB.X.Y (default: off)",
    ),
    parse_FL_coverage: bool = typer.Option(
        False,
        help="Parse `full_length_coverage=` from sequence ID.",
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    parse_matchAnnot(fa_or_fq, match_filename, not_pbid, parse_FL_coverage)


if __name__ == "__main__":
    typer.run(main)
