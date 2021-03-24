#!/usr/bin/env python
import sys
from collections import defaultdict

from Bio import SeqIO
from cupcake.logging import cupcake_logger as logger


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


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser("Parse MatchAnnot result")
    parser.add_argument(
        "fa_or_fq",
        help="Fasta/Fastq filename used to create the SAM file for matchAnnot",
    )
    parser.add_argument("match_filename", help="MatchAnnot filename")
    parser.add_argument(
        "--not_pbid",
        action="store_true",
        default=False,
        help="Turn this on if not sequence ID is not PB.X.Y (default: off)",
    )
    parser.add_argument(
        "--parse_FL_coverage",
        action="store_true",
        default=False,
        help="Parse `full_length_coverage=` from sequence ID.",
    )

    args = parser.parse_args()
    parse_matchAnnot(
        args.fa_or_fq, args.match_filename, args.not_pbid, args.parse_FL_coverage
    )


if __name__ == "__main__":
    main()
