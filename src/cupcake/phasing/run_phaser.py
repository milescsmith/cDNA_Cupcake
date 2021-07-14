#!/usr/bin/env python

__author__ = "etseng@pacb.com"

import sys
from enum import Enum
from pathlib import Path
from subprocess import check_output

import typer
from Bio import SeqIO

from cupcake import version_callback
from cupcake import cupcake_logger as logger
from cupcake.phasing.io import VariantPhaseCleaner, VariantPhaser
from cupcake.phasing.io.MPileUpVariantCaller import MPileUPVariant
from cupcake.phasing.io.SAMMPileUpReader import MPileUpReader

# try:
#     import vcf
# except ImportError:
#     print("Cannot import vcf! Please install pyvcf!", file=sys.stderr)
#     sys.exit(-1)

MIN_COVERAGE = 10  # minimum number of FL reads for a gene to do SNP calling and phasing
ERR_SUB = 0.005
MAX_DIFF_ALLOWED = 3  # maximum difference in bases allowed for two haplotype strings
MIN_PERC_ALLOWED = 0.25  # minimum percent of total count for an allele, can be adjusted by ploidy (ex: n=6, means this goes down to 1/6)
PVAL_CUTOFF = 0.01
MIN_AF_AT_ENDS = 0.10  # minimum minor allele frequency for SNPs at ends, which tend to have unreliable alignments


class strand_direction(str, Enum):
    plus = "+"
    minus = "-"


app = typer.Typer(name="cupcake.phasing.run_phase")


def set_to_kill(
    fastx_filename: str,
    sam_filename: str,
    mpileup_filename: str,
    read_stat: str,
    mapping_filename: str,
    output_prefix: str,
    strand: str,
    partial_ok: bool,
    pval_cutoff: float,
    ploidy: int,
) -> None:
    # remove potential past run output
    past_files = [
        f"{output_prefix}.NO_SNPS_FOUND",
        f"{output_prefix}.NO_HAPS_FOUND",
        f"{output_prefix}.log",
        f"{output_prefix}.human_readable.txt",
        f"{output_prefix}.vcf",
        f"{output_prefix}.cleaned.human_readable.txt",
        f"{output_prefix}.cleaned.vcf",
    ]

    for file in past_files:
        f = Path(file)
        if f.exists():
            f.unlink()

    # (1) read the mpileup and vall variants
    reader = MPileUpReader(mpileup_filename)
    recs = list(reader)
    vc = MPileUPVariant(
        recs,
        min_cov=MIN_COVERAGE,
        err_sub=ERR_SUB,
        expected_strand=strand,
        pval_cutoff=pval_cutoff,
    )
    vc.call_variant()
    print(vc.variant)

    if len(vc.variant) == 0:
        check_output(["touch", f"{output_prefix}.NO_SNPS_FOUND"])
        logger.critical("No SNPs found. END.")
        sys.exit(0)

    # (2) for each CCS read, assign a haplotype (or discard if outlier)
    pp = VariantPhaser.VariantPhaser(vc)
    pp.phase_variant(sam_filename, fastx_filename, output_prefix, partial_ok=partial_ok)
    pp.haplotypes
    pp.haplotypes.get_haplotype_vcf_assignment()

    # (3) phase isoforms
    seqids = {
        r.id
        for r in SeqIO.parse(
            open(fastx_filename), VariantPhaser.type_fa_or_fq(fastx_filename)
        )
    }
    isoform_tally = VariantPhaser.phase_isoforms(read_stat, seqids, pp)
    if len(isoform_tally) == 0:
        check_output(["touch", f"{output_prefix}.NO_HAPS_FOUND"])
        logger.critical("No good haps found. END.")
        sys.exit(0)
    pp.haplotypes.write_haplotype_to_vcf(mapping_filename, isoform_tally, output_prefix)

    # (4) clean isoforms
    hap_count = VariantPhaseCleaner.make_haplotype_counts(isoform_tally)

    # (5) error correct haplotypes
    #  if diploid, use exhaustive search
    #  otherwise, use hap counts (ToDo: make this work with exhaustive search later)
    variants = [
        [base.upper() for base, count in vc.variant[pos]] for pos in pp.accepted_pos
    ]

    if ploidy == 2 and all(len(_) == 2 for _ in variants):
        (
            diff_arr,
            hap_count_ordered,
        ) = VariantPhaseCleaner.infer_haplotypes_via_exhaustive_diploid_only(
            pp.haplotypes, variants
        )
    else:
        min_perc_allowed = min(
            MIN_PERC_ALLOWED, 1 / (ploidy + 4)
        )  # this is a heuristic, allowing for some alleles to be much less expressde than others
        diff_arr, hap_count_ordered = VariantPhaseCleaner.infer_haplotypes_via_min_diff(
            pp.haplotypes.haplotypes,
            hap_count,
            ploidy,
            MAX_DIFF_ALLOWED,
            min_perc_allowed,
        )

    if diff_arr is None:
        check_output(["touch", f"{output_prefix}.cleaned.NO_HAPS_FOUND"])
        logger.critical("No good haps found. END.")
        sys.exit(0)

    m, new_hap, new_isoform_tally = VariantPhaseCleaner.error_correct_haplotypes(
        pp.haplotypes, isoform_tally, diff_arr, hap_count_ordered
    )
    # write out the mapping relationship between: FL CCS --> (pre-corrected) hap --> error-corrected hap
    with open(f"{output_prefix}.cleaned.hap_info.txt", "w") as f:
        f.write("id,hap_preclean,hap_postclean\n")
        for seqid, old_i in pp.seq_hap_info.items():
            f.write(
                f"{seqid},{pp.haplotypes.haplotypes[old_i]},{new_hap.haplotypes[m[old_i]]}\n"
            )

    new_hap.get_haplotype_vcf_assignment()
    new_hap.write_haplotype_to_vcf(
        mapping_filename, new_isoform_tally, f"{output_prefix}.cleaned"
    )


@app.command(name="")
def main(
    fastx_filename: str = typer.Argument(..., help="Input FLNC fasta or fastq"),
    sam_filename: str = typer.Argument(...),
    mpileup_filename: str = typer.Argument(
        ...,
    ),
    read_stat: str = typer.Argument(
        ...,
    ),
    mapping_filename: str = typer.Argument(...),
    output_prefix: str = typer.Option(..., "--output_prefix", "-o"),
    strand: strand_direction = typer.Option(...),
    partial_ok: bool = typer.Option(
        False,
    ),
    pval_cutoff: float = typer.Option(
        PVAL_CUTOFF,
        "--pval_cutoff",
        "-p",
    ),
    ploidy: int = typer.Option(2),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
    # parser.add_argument("-e", "--err_sub", default=ERR_SUB, type=float, help="Estimated substitution error rate (default: 0.005)")
) -> None:
    set_to_kill(
        fastx_filename=fastx_filename,
        sam_filename=sam_filename,
        mpileup_filename=mpileup_filename,
        read_stat=read_stat,
        mapping_filename=mapping_filename,
        output_prefix=output_prefix,
        strand=strand,
        partial_ok=partial_ok,
        pval_cutoff=pval_cutoff,
        ploidy=ploidy,
    )


if __name__ == "main":
    typer.run(main)
