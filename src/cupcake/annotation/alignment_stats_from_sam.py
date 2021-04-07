#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Given an input fasta/fastq and alignment SAM (from minimap2, GMAP, STAR, etc),
Output two report files.

Report #1. For each alignment:
  (1) sequence id
  (2) alignment coverage (just for this alignment)
  (3) alignmen identity
  (4) number of mismmatches
  (5) number of insertions
  (6) number of deletions
  (7) number of exons

Report #2. For each junction:
  (1) sequence id
  (2) donor position (<chr>:<strand>:<0-based pos>)
  (3) donor sequence ("GT")
  (4) closest known donor site (NA if no ref, 0 means right on spot)
  (5) acceptor position (<chr>:<strand>:<0-based pos>)
  (6) acceptor sequence ("AG")
  (7) closest known donor site (NA if no ref)
"""
import bisect
from collections import defaultdict
from csv import DictWriter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import typer
from Bio import SeqIO
from Bio.Seq import Seq
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence.BioReaders import GMAPSAMReader
from cupcake.sequence.GFF import collapseGFFReader

FIELDNAMES_REPORT1 = [
    "seqid",
    "coverage",
    "identity",
    "num_sub",
    "num_ins",
    "num_del",
    "num_exons",
]
FIELDNAMES_REPORT2 = [
    "seqid",
    "donor_pos",
    "donor_seq",
    "donor_dist",
    "acceptor_pos",
    "acceptor_seq",
    "acceptor_dist",
]


app = typer.Typer(name="alignment_stats_from_sam")


def type_fa_or_fq(file: str) -> str:
    file = file.upper()
    if file.endswith(".FA") or file.endswith(".FASTA"):
        return "fasta"
    else:
        return "fastq"


def read_annotation_for_junction_info(gff_filename: str) -> defaultdict:
    """
    :param gff_filename: annotation GFF
    :return: dict of (chrom, strand, 'donor' or 'acceptor') --> sorted list of donor or acceptor site. all 0-based.
    """
    d = defaultdict(lambda: set())
    for r in collapseGFFReader(gff_filename):
        if r.strand == "+":
            for i in range(0, len(r.ref_exons) - 1):
                d[(r.chr, r.strand, "donor")].add(r.ref_exons[i].end - 1)
                d[(r.chr, r.strand, "acceptor")].add(r.ref_exons[i + 1].start)
        else:
            for i in range(0, len(r.ref_exons) - 1):
                d[(r.chr, r.strand, "acceptor")].add(r.ref_exons[i].end - 1)
                d[(r.chr, r.strand, "donor")].add(r.ref_exons[i + 1].start)
    for k in d:
        d[k] = list(d[k])
        d[k].sort()
    return d


def get_closest_junction_dist(
    junction_info: Dict[Tuple[str, str, str], List[int]],
    chrom: str,
    strand: str,
    pos1: int,
    pos2: int,
) -> Tuple[int, int]:
    """
    pos1, pos2 would be the donor/acceptor, in 0-based index
    always pos1 < pos2
    if strand == '+': pos1 is the donor
    if strand == '-': pos1 is the acceptor

    :return: (dist to closest donor, dist to closest acceptor)
    """

    def get_min_dist(distlist: List[int], index: int, pos: int) -> int:
        if index == 0:
            return distlist[0] - pos
        elif index == len(distlist):
            return pos - distlist[-1]
        else:
            return min(distlist[index] - pos, pos - distlist[index - 1])

    if strand == "+":
        distlist1 = junction_info[(chrom, strand, "donor")]
        distlist2 = junction_info[(chrom, strand, "acceptor")]
        if (
            len(distlist1) == 0 or len(distlist2) == 0
        ):  # the list is empty, no junctions
            return "NA", "NA"
        else:
            i = bisect.bisect_left(distlist1, pos1)
            j = bisect.bisect_left(distlist2, pos2)
            return get_min_dist(distlist1, i, pos1), get_min_dist(distlist2, j, pos2)
    else:
        distlist1 = junction_info[(chrom, strand, "acceptor")]
        distlist2 = junction_info[(chrom, strand, "donor")]
        if (
            len(distlist1) == 0 or len(distlist2) == 0
        ):  # the list is empty, no junctions
            return "NA", "NA"
        else:
            i = bisect.bisect_left(distlist1, pos1)
            j = bisect.bisect_left(distlist2, pos2)
            return get_min_dist(distlist2, j, pos2), get_min_dist(distlist1, i, pos1)


def get_donor_acceptor(
    genome_d: Dict[str, Seq], chrom: str, strand: str, pos1: int, pos2: int
) -> Tuple[str, str]:
    """
    pos1, pos2 would be the donor/acceptor, in 0-based index
    always pos1 < pos2
    if strand == '+': pos1 is the donor
    if strand == '-': pos1 is the acceptor

    :return: donor sequence (ex: "GT"), acceptor sequence (ex: "AG")
    """
    seq1 = genome_d[chrom][(pos1 + 1) : (pos1 + 3)].seq
    seq2 = genome_d[chrom][(pos2 - 2) : pos2].seq
    if strand == "+":
        return str(seq1).upper(), str(seq2).upper()
    else:
        return (
            str(seq2.reverse_complement()).upper(),
            str(seq1.reverse_complement()).upper(),
        )


def evaluate_alignment_sam(
    input_fa_or_fq: str,
    sam_filename: str,
    genome_d: Dict[str, Seq],
    output_prefix: str,
    junction_info=None,
) -> None:
    h1 = Path(f"{output_prefix}.alignment_report.txt")
    h2 = Path(f"{output_prefix}.junction_report.txt")

    w1 = DictWriter(open(h1, "w"), fieldnames=FIELDNAMES_REPORT1)
    w2 = DictWriter(open(h2, "w"), fieldnames=FIELDNAMES_REPORT2)
    w1.writeheader()
    w2.writeheader()

    # FIELDNAMES_REPORT1 = ['seqid', 'coverage', 'identity', 'num_sub', 'num_ins', 'num_del', 'num_exons']
    # FIELDNAMES_REPORT2 = ['seqid', 'donor_pos', 'donor_seq', 'donor_dist', 'acceptor_pos', 'acceptor_seq', 'acceptor_dist']

    query_len_dict = {
        r.id: len(r.seq)
        for r in SeqIO.parse(open(input_fa_or_fq), type_fa_or_fq(input_fa_or_fq))
    }
    for r in GMAPSAMReader(sam_filename, True, query_len_dict=query_len_dict):
        if r.sID == "*":  # unaligned
            rec1 = {
                "seqid": r.qID,
                "coverage": "NA",
                "identity": "NA",
                "num_sub": "NA",
                "num_ins": "NA",
                "num_del": "NA",
                "num_exons": "NA",
            }
            w1.writerow(rec1)
            continue
        rec1 = {
            "seqid": r.qID,
            "coverage": r.qCoverage,
            "identity": r.identity,
            "num_sub": r.num_nonmatches - r.num_del - r.num_ins,
            "num_ins": r.num_ins,
            "num_del": r.num_del,
            "num_exons": len(r.segments),
        }
        w1.writerow(rec1)
        for i in range(0, len(r.segments) - 1):
            rec2 = {"seqid": r.qID}
            seq1, seq2 = get_donor_acceptor(
                genome_d,
                r.sID,
                r.flag.strand,
                r.segments[i].end - 1,
                r.segments[i + 1].start,
            )
            if r.flag.strand == "+":
                rec2["donor_pos"] = f"{r.sID}:+:{r.segments[i].end - 1}"
                rec2["acceptor_pos"] = f"{r.sID}:+:{r.segments[i + 1].start}"
            else:
                rec2["donor_pos"] = f"{r.sID}:-:{r.segments[i + 1].start}"
                rec2["acceptor_pos"] = f"{r.sID}:-:{r.segments[i].end - 1}"
            rec2["donor_seq"] = seq1
            rec2["acceptor_seq"] = seq2
            if junction_info is not None:
                rec2["donor_dist"], rec2["acceptor_dist"] = get_closest_junction_dist(
                    junction_info,
                    r.sID,
                    r.flag.strand,
                    r.segments[i].end - 1,
                    r.segments[i + 1].start,
                )
            else:
                rec2["donor_dist"] = "NA"
                rec2["acceptor_dist"] = "NA"
            w2.writerow(rec2)


@app.command(name="")
def main(
    input_file: str = typer.Option(..., "--input", "-i", help="Input fasta or fastq."),
    sam_filename: str = typer.Option(
        ..., "--sam_filename", "-s", help="Aligned SAM filename."
    ),
    genome_filename: str = typer.Option(
        ..., "--genome_filename", "-g", help="Genome fasta."
    ),
    output_prefix: str = typer.Option(
        ..., "--output_prefix", "-o", help="Output prefix."
    ),
    gff: Optional[str] = typer.Option(None, "--gff", help="Annotation GFF."),
):

    # read genome
    logger.info(f"Reading genome {genome_filename}...")
    genome_d = SeqIO.to_dict(SeqIO.parse(open(genome_filename), "fasta"))

    # read gff
    if gff is not None:
        logger.info(f"Reading annotation {gff}...")
        junction_info = read_annotation_for_junction_info(gff)
    else:
        junction_info = None

    evaluate_alignment_sam(
        input_file, sam_filename, genome_d, output_prefix, junction_info
    )


if __name__ == "__main__":
    typer.run(main)
