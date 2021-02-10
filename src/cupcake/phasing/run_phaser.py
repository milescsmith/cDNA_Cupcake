#!/usr/bin/env python

__author__ = "etseng@pacb.com"

import sys
from argparse import ArgumentParser
from subprocess import check_output
from pathlib import Path

from Bio import SeqIO

from cupcake.phasing.io.MPileUpVariantCaller import MPileUPVariant
from cupcake.phasing.io.SAMMPileUpReader import MPileUpReader
from cupcake.phasing.io import VariantPhaseCleaner, VariantPhaser

# try:
#     import vcf
# except ImportError:
#     print("Cannot import vcf! Please install pyvcf!", file=sys.stderr)
#     sys.exit(-1)

MIN_COVERAGE = 10  # minimum number of FL reads for a gene to do SNP calling and phasing
ERR_SUB = 0.005
MAX_DIFF_ALLOWED = 3  # maximum difference in bases allowed for two haplotype strings
MIN_PERC_ALLOWED = (
    0.25
)  # minimum percent of total count for an allele, can be adjusted by ploidy (ex: n=6, means this goes down to 1/6)
PVAL_CUTOFF = 0.01
MIN_AF_AT_ENDS = (
    0.10
)  # minimum minor allele frequency for SNPs at ends, which tend to have unreliable alignments


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
    recs = [r for r in reader]
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
        print("No SNPs found. END.", file=sys.stderr)
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
        print("No good haps found. END.", file=sys.stderr)
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

    if ploidy == 2 and all(len(vars) == 2 for vars in variants):
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
        print("No good haps found. END.", file=sys.stderr)
        sys.exit(0)

    m, new_hap, new_isoform_tally = VariantPhaseCleaner.error_correct_haplotypes(
        pp.haplotypes, isoform_tally, diff_arr, hap_count_ordered
    )
    # write out the mapping relationship between: FL CCS --> (pre-corrected) hap --> error-corrected hap
    with open(output_prefix + ".cleaned.hap_info.txt", "w") as f:
        f.write("id,hap_preclean,hap_postclean\n")
        for seqid, old_i in pp.seq_hap_info.items():
            f.write(
                f"{seqid},{pp.haplotypes.haplotypes[old_i]},{new_hap.haplotypes[m[old_i]]}\n"
            )

    new_hap.get_haplotype_vcf_assignment()
    new_hap.write_haplotype_to_vcf(
        mapping_filename, new_isoform_tally, f"{output_prefix}.cleaned"
    )


def main():
    parser = ArgumentParser()
    parser.add_argument("fastx_filename", type=str, help="Input FLNC fasta or fastq")
    parser.add_argument("sam_filename", type=str)
    parser.add_argument("mpileup_filename", type=str)
    parser.add_argument("read_stat", type=str)
    parser.add_argument("mapping_filename", type=str)
    parser.add_argument("-o", "--output_prefix", type=str, required=True)
    parser.add_argument("--strand", choices=["+", "-"], type=str, required=True)
    parser.add_argument("--partial_ok", default=False, type=bool, action="store_true")
    parser.add_argument("-p", "--pval_cutoff", default=PVAL_CUTOFF, type=float)
    parser.add_argument("-n", "--ploidy", type=int, default=2)
    # parser.add_argument("-e", "--err_sub", default=ERR_SUB, type=float, help="Estimated substitution error rate (default: 0.005)")

    args = parser.parse_args()
    set_to_kill(
        fastx_filename=args.fastx_filename,
        sam_filename=args.sam_filename,
        mpileup_filename=args.mpileup_filename,
        read_stat=args.read_stat,
        mapping_filename=args.mapping_filename,
        output_prefix=args.output_prefix,
        strand=args.strand,
        partial_ok=args.partial_ok,
        pval_cutoff=args.pval_cutoff,
        ploidy=args.ploidy,
    )


if __name__ == "main":
    main()
