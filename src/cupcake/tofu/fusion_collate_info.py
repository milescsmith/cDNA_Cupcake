#!/usr/bin/env python

import re
import sys
from collections import defaultdict
from csv import DictReader, DictWriter
from pathlib import Path
from typing import Optional, Tuple

import typer
from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger
from cupcake.sequence.GFF import collapseGFFReader

fusion_pbid = re.compile(r"PBfusion.(\d+).(\d+)")
"""
Run after fusion_finder.py + SQANTI3 classification
"""

FIELDS = [
    "UniqueID",
    "FusionName",
    "LeftGeneName",
    "LeftGeneID",
    "LeftBreakpoint",
    "LeftFlankingSequence",
    "RightGeneName",
    "RightGeneID",
    "RightBreakpoint",
    "RightFlankingSequence",
    "JunctionSupport",
    "SpanningReads",
    "ReadCountScore",
    "Sequence",
    "LeftORF",
    "RightORF",
    "LeftExonCount",
    "RightExonCount",
    "LeftCDSExonCount",
    "RightCDSExonCount",
    "Comments",
]


app = typer.Typer(name="cupcake.tofu.fusion_collate_into")


def get_breakpoint_n_seq(
    r1: dict, r2: dict, genome_dict: Optional[str] = None, flanking_size: int = 50
) -> Tuple[str, str, str, str]:
    if r1.strand == "+":
        left_breakpoint = f"{r1.chr}:{r1.end}:+"
        if genome_dict is not None:
            left_seq = str(genome_dict[r1.chr][r1.end - flanking_size : r1.end].seq)
        else:
            left_seq = "NA"
    else:
        left_breakpoint = f"{r1.chr}:{r1.start+1}:-"
        if genome_dict is not None:
            left_seq = str(
                genome_dict[r1.chr][r1.start : r1.start + flanking_size]
                .reverse_complement()
                .seq
            )
        else:
            left_seq = "NA"
    if r2.strand == "+":
        right_breakpoint = f"{r2.chr}:{r2.start}:+"
        if genome_dict is not None:
            right_seq = str(
                genome_dict[r2.chr][r2.start : r2.start + flanking_size].seq
            )
        else:
            right_seq = "NA"
    else:
        right_breakpoint = f"{r2.chr}:{r2.end}:-"
        if genome_dict is not None:
            right_seq = str(
                genome_dict[r2.chr][r2.end - flanking_size : r2.end]
                .reverse_complement()
                .seq
            )
        else:
            right_seq = "NA"
    return left_breakpoint, left_seq, right_breakpoint, right_seq


def collate_info(
    fusion_prefix: str,
    class_filename: str,
    genepred_filename: str,
    total_fl_count: Optional[int] = None,
    config_filename: Optional[str] = None,
    genome_dict: Optional[dict] = None,
    cds_gff_filename: Optional[str] = None,
    min_fl_count: int = 2,
    min_breakpoint_dist_kb: int = 10,
    include_Mt_genes: bool = False,
) -> None:

    global_info = {}  # holding information for general information
    if config_filename is not None:
        logger.info(f"Reading config file {config_filename}...")
        for line in open(config_filename):
            k, v = line.strip().split("=")
            global_info[k] = v

    gene_to_id = {}  # gene name --> ensembl ID
    for line in open(genepred_filename):
        raw = line.strip().split()
        gene_to_id[raw[11]] = raw[0]

    d = defaultdict(lambda: {})  # PBfusion.X --> isoform index -> sqanti3 record
    orf_dict = {}
    # read SQANTI3 classification file
    for r in DictReader(open(class_filename), delimiter="\t"):
        m = fusion_pbid.match(r["isoform"])
        if m is None:
            logger.error("ERROR: fusion pbid must follow format `PBfusion.X.Y`. Abort!")
            sys.exit(-1)
        gene_index, isoform_index = m.group(1), m.group(2)
        d[gene_index][isoform_index] = r
        orf_dict[r["isoform"]] = r["ORF_seq"]

    # get sequences
    seq_dict = {
        r.id.split("|")[0]: r.seq
        for r in SeqIO.parse(open(f"{fusion_prefix}.rep.fa"), "fasta")
    }

    # get count information
    count_d = defaultdict(lambda: "NA")
    count_filename = f"{fusion_prefix}.abundance.txt"
    if Path(count_filename).exists():
        for r in DictReader(open(count_filename), delimiter="\t"):
            count_d[r["pbid"]] = int(r["count_fl"])

    if total_fl_count is None:
        logger.info(
            "Total FL count not given --- using the sum FL count from fusions only instead."
        )
        total_fl_count = sum(count_d.values())

    # get breakpoint information
    gff_d = defaultdict(lambda: {})  # PBfusion.X --> isoform index -> sqanti3 record
    if cds_gff_filename is None:
        gff_filename = f"{fusion_prefix}.gff"
    else:
        gff_filename = cds_gff_filename

    for r in collapseGFFReader(gff_filename):
        m = fusion_pbid.match(r.seqid)
        if m is None:
            logger.error(
                f"ERROR: fusion pbid in {gff_filename} must follow format `PBfusion.X.Y`. Abort!"
            )
            sys.exit(-1)
        gene_index, isoform_index = m.group(1), int(m.group(2))
        gff_d[gene_index][isoform_index] = r
        if r.strand not in ("+", "-"):
            logger.error(
                f"ERROR: fusion {r.seqid} did not specify strand in {gff_filename}! Abort!"
            )
            sys.exit(-1)

    fields2 = list(global_info.keys()) + FIELDS
    with open(f"{fusion_prefix}.annotated.txt", "w") as f, open(
        f"{fusion_prefix}.annotated_ignored.txt", "w"
    ) as f_bad:
        writer = DictWriter(f, fields2, delimiter=",")
        writer.writeheader()
        writer_bad = DictWriter(f_bad, fields2, delimiter=",")
        writer_bad.writeheader()

        for gene_index, iso_dict in d.items():
            iso_dict = list(iso_dict.items())  # (isoform index, classification record)
            iso_dict.sort(key=lambda x: x[0])
            has_novel = any(
                r["associated_gene"].startswith("novelGene")
                or r["associated_gene"] == ""
                for junk, r in iso_dict
            )
            pbid = f"PBfusion.{str(gene_index)}"

            gff_info = list(gff_d[gene_index].items())
            gff_info.sort(key=lambda x: x[0])

            rec1 = gff_info[0][1]
            rec2 = gff_info[-1][1]
            (
                left_breakpoint,
                left_seq,
                right_breakpoint,
                right_seq,
            ) = get_breakpoint_n_seq(rec1, rec2, genome_dict)
            left_exon_count = len(rec1.ref_exons)
            right_exon_count = len(rec2.ref_exons)
            gene1 = iso_dict[0][1]["associated_gene"]
            gene2 = iso_dict[-1][1]["associated_gene"]

            if cds_gff_filename is not None:
                left_cds_exon_count = len(rec1.cds_exons)
                right_cds_exon_count = len(rec2.cds_exons)
            else:
                left_cds_exon_count = "NA"
                right_cds_exon_count = "NA"

            left_orf, right_orf = "NA", "NA"
            if orf_dict is not None:
                seqid1 = gff_info[0][1].seqid
                seqid2 = gff_info[-1][1].seqid
                left_orf = orf_dict[seqid1]
                right_orf = orf_dict[seqid2]

            info = {
                "UniqueID": pbid,
                "FusionName": "--".join(
                    [_r["associated_gene"] for (_index, _r) in iso_dict]
                ),
                "LeftGeneName": gene1,
                "LeftGeneID": gene_to_id[gene1] if gene1 in gene_to_id else "NA",
                "LeftBreakpoint": left_breakpoint,
                "LeftFlankingSequence": left_seq,
                "RightGeneName": gene2,
                "RightGeneID": gene_to_id[gene2] if gene2 in gene_to_id else "NA",
                "RightBreakpoint": right_breakpoint,
                "RightFlankingSequence": right_seq,
                "JunctionSupport": "NA",
                "SpanningReads": count_d[pbid],
                "ReadCountScore": (count_d[pbid] * (10 ** 6) / total_fl_count)
                if count_d[pbid] != "NA"
                else "NA",
                "Sequence": seq_dict[pbid],
                "LeftORF": left_orf,
                "RightORF": right_orf,
                "LeftExonCount": left_exon_count,
                "RightExonCount": right_exon_count,
                "LeftCDSExonCount": left_cds_exon_count,
                "RightCDSExonCount": right_cds_exon_count,
                "Comments": "PASS",
            }
            info.update(global_info)

            left_chr, left_break, left_strand = left_breakpoint.split(":")
            right_chr, right_break, right_strand = right_breakpoint.split(":")

            if has_novel:
                info["Comments"] = "FAIL:NovelGene"
            elif gene1 == gene2:
                info["Comments"] = "FAIL:SameGene"
            elif info["SpanningReads"] != "NA" and info["SpanningReads"] < min_fl_count:
                info["Comments"] = "FAIL:TooFewFLReads"
            elif not include_Mt_genes and (
                gene1.startswith("MT-") or gene2.startswith("MT-")
            ):
                info["Comments"] = "FAIL:MtGenes"
            elif (
                left_chr == right_chr
                and abs(int(left_break) - int(right_break)) / 1000
                <= min_breakpoint_dist_kb
            ):
                info["Comments"] = "FAIL:BreakpointTooClose"

            if info["Comments"].startswith("FAIL:"):
                writer_bad.writerow(info)
            else:
                writer.writerow(info)


@app.command(name="")
def main(
    fusion_prefix: str = typer.Argument(
        ..., help="Prefix for fusion files (ex: my.fusion)"
    ),
    class_filename: str = typer.Argument(..., help="SQANTI3 classification filename"),
    genepred_filename: str = typer.Argument(
        ..., help="GenePred filename used by SQANTI3 classification"
    ),
    cds_gff: Optional[str] = typer.Option(None, help="CDS GFF filename"),
    total_fl_count: Optional[int] = typer.Option(
        None, help="Total FL count used to normalize fusion counts"
    ),
    config: Optional[str] = typer.Option(
        None, help="(optional) Additional information to include in the output"
    ),
    genome: Optional[str] = typer.Option(None, help="Reference genome"),
    min_fl_count: int = typer.Option(2, help="Minimum FL count"),
    min_breakpoint_dist_kb: int = typer.Option(
        10, help="Minimum breakpoint distance, in kb)"
    ),
    include_Mt_genes: bool = typer.Option(False, help="Include Mt genes"),
) -> None:
    if genome is not None:
        genome_dict = SeqIO.to_dict(SeqIO.parse(open(genome), "fasta"))
        print(f"Finished reading reference genome {genome}.")
    else:
        genome_dict = None

    collate_info(
        fusion_prefix,
        class_filename,
        genepred_filename,
        total_fl_count=total_fl_count,
        config_filename=config,
        genome_dict=genome_dict,
        cds_gff_filename=cds_gff,
        min_fl_count=min_fl_count,
        min_breakpoint_dist_kb=min_breakpoint_dist_kb,
        include_Mt_genes=include_Mt_genes,
    )


if __name__ == "__main__":
    typer.run(main)
