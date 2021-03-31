__author__ = "etseng@pacb.com"
import glob
from collections import defaultdict
from csv import DictReader, DictWriter
from pathlib import Path

import typer
import vcfpy
from bx.intervals import IntervalTree
from cupcake.logging import cupcake_logger as logger
from cupcake.phasing.io.SAMMPileUpReader import MPileUpReader

app = typer.Typer(name="cupcake.phasing.util.evaluate_snp_with_genome")


def read_fake_mapping(fake_genome_mapping_filename):
    """
    :param fake_genome_mapping_filename: fake genome mapping file
    :return: dict of fake position -> (chr, ref position)
    """
    fake_map = {}
    with open(fake_genome_mapping_filename) as (f):
        for line in f:
            fake_pos, ref_chr, ref_pos = line.strip().split(",")
            fake_map[int(fake_pos)] = (ref_chr, int(ref_pos))

    return fake_map


def get_positions_to_recover(
    fake_genome_mapping_filename, mpileup_filename, genome_snp, min_cov=40
):
    """
    Return a list of genomic positions that isophase should call SNPs for. Must be:

    1. covered (based on ccs.mpileup) with at least MIN_COVERAGE (default: 40)
    2. is a SNP (based on genome VCF)
    """
    fake_map = read_fake_mapping(fake_genome_mapping_filename)
    good_positions = []
    cov_at_pos = defaultdict(
        lambda: 0
    )  # dict of (chrom,pos) --> (coverage)  // this is used for reports later
    for rec in MPileUpReader(mpileup_filename):
        if rec is None:  # means coverage zero, ignore
            continue
        genome_chr, genome_pos = fake_map[rec.pos]
        cov_at_pos[(genome_chr, genome_pos)] = rec.nCov
        if rec.nCov >= min_cov:
            genome_chr = genome_chr.split("|")[0]
            if (
                genome_chr in genome_snp
                and genome_pos in genome_snp[genome_chr]
                and genome_snp[genome_chr][genome_pos].is_snp
            ):
                good_positions.append((genome_chr, genome_pos))

    return good_positions, cov_at_pos


def eval_isophase(
    isophase_vcf,
    genome_snp,
    good_positions,
    cov_at_pos,
    repeat_by_chrom,
    shortread_cov,
    writer_f,
    name="NA",
    strand="NA",
):
    for r in vcfpy.Reader(isophase_vcf):
        out = {
            "dir": name,
            "chrom": "NA",
            "pos": r.POS,
            "strand": strand,
            "ref": r.REF,
            "alt_Short": "NA",
            "alt_PB": "NA",
            "in_Short": "NA",
            "in_PB": "NA",
            "cov_Short": "NA",
            "cov_PB": "NA",
            "genomic_HP": "NA",
        }

        r.CHROM = r.CHROM.split("|")[0]
        out["chrom"] = r.CHROM
        out["alt_PB"] = r.ALT[0]

        out["genomic_HP"] = (
            "Y"
            if (
                r.CHROM in repeat_by_chrom
                and len(repeat_by_chrom[r.CHROM].find(r.POS, r.POS)) > 0
            )
            else "N"
        )
        try:
            out["cov_Short"] = shortread_cov[r.CHROM][r.POS]
        except KeyError:
            out["cov_Short"] = 0
        out["cov_PB"] = cov_at_pos[r.CHROM, r.POS - 1]
        if (r.CHROM, r.POS) not in good_positions:
            out["alt_Short"] = "NA"
            out["in_Short"] = "N"
            out["in_PB"] = "Y"
        else:
            out["alt_Short"] = genome_snp[r.CHROM][r.POS].ALT[0]
            out["in_Short"] = "Y"
            out["in_PB"] = "Y"
            good_positions.remove((r.CHROM, r.POS))
        writer_f.writerow(out)

    # now we write out everything that is only in Shortread
    for chrom, pos in good_positions:
        out = {
            "dir": name,
            "chrom": chrom,
            "pos": pos,
            "strand": strand,
            "ref": genome_snp[chrom][pos].REF,
            "alt_Short": genome_snp[chrom][pos].ALT[0],
            "alt_PB": "NA",
            "in_Short": "Y",
            "in_PB": "N",
            "cov_Short": "NA",
            "cov_PB": cov_at_pos[chrom, pos - 1],
            "genomic_HP": "Y"
            if (
                chrom in repeat_by_chrom
                and len(repeat_by_chrom[chrom].find(pos, pos)) > 0
            )
            else "N",
        }
        try:
            out["cov_Short"] = shortread_cov[chrom][pos]
        except KeyError:
            out["cov_Short"] = 0
        writer_f.writerow(out)


def brangus(vcf_filename, out_filename, unzip_snps=None):
    if unzip_snps is None:
        unzip_snps = defaultdict(lambda: {})
        for r in vcfpy.Reader(vcf_filename):
            unzip_snps[r.CHROM][r.POS] = r

    logger.info(f"Finished reading {vcf_filename}")
    with open(out_filename, "w") as out_f:
        FIELDS = [
            "dir",
            "chrom",
            "pos",
            "strand",
            "ref",
            "alt_Short",
            "alt_PB",
            "in_Short",
            "in_PB",
            "cov_Short",
            "cov_PB",
            "genomic_HP",
        ]
        writer = DictWriter(out_f, FIELDS, delimiter="\t")
        writer.writeheader()
        dirs = glob.glob("by_loci/*size*/")
        for d1 in dirs:
            mpileup = Path(d1, "ccs.mpileup")
            mapfile = Path(d1, "fake.mapping.txt")
            vcffile = Path(d1, "phased.partial.vcf")
            config = Path(d1, "config")
            nosnp = Path(d1, "phased.partial.NO_SNPS_FOUND")
            if not vcffile.exists():
                if not nosnp.exists():
                    logger.error(f"Skipping {d1} because no SNPs found.")
            else:
                logger.info(f"Evaluating {d1}.")
                strand = "NA"
                if config.exists():  # find the strand this gene family is on
                    for line in open(config):
                        if line.startswith("ref_strand="):
                            strand = line.strip().split("=")[1]
                good_positions, cov_at_pos = get_positions_to_recover(
                    mapfile, mpileup, unzip_snps, min_cov=30
                )
                name = d1.split("/")[1]
                eval_isophase(
                    vcffile,
                    unzip_snps,
                    good_positions,
                    cov_at_pos,
                    {},
                    {},
                    writer,
                    name,
                    strand,
                )

    return


def main_maize(ki11_snps=None, dirs=None):
    if ki11_snps is None:
        ki11_snps = defaultdict(lambda: {})  # chrom -> pos -> VCF record
        debug_count = 0
        for r in vcfpy.Reader("B73Ki11.q20.vcf"):
            ki11_snps[r.CHROM][r.POS] = r
            # if debug_count > 100000: break
            debug_count += 1

    logger.info("Finished reading B73Ki11.q20.vcf.")

    ki11_shortread_cov = defaultdict(lambda: {})  # chrom -> pos -> short read cov
    # read the raw Ki11 pileup to get coverage in places where no SNPs were called
    for r in MPileUpReader("Ki11.raw.mpileup"):
        if r is not None:
            ki11_shortread_cov[r.chr][r.pos] = r.cov
    logger.info("Fnished reading Ki11.raw.mpileup.")

    repeat_by_chrom = {}
    # read the Tandem Repeat Finder summary
    for r in DictReader(open("B73_RefV4.fa.repeat_list.txt"), delimiter="\t"):
        if r["chrom"] not in repeat_by_chrom:
            repeat_by_chrom[r["chrom"]] = IntervalTree()
        repeat_by_chrom[r["chrom"]].add(int(r["start0"]), int(r["end1"]))

    logger.info("Finished reading B73_RefV4.fa.repeat_list.txt.")

    FIELDS = [
        "dir",
        "chrom",
        "pos",
        "ref",
        "alt_Short",
        "alt_PB",
        "in_Short",
        "in_PB",
        "cov_Short",
        "cov_PB",
        "genomic_HP",
    ]
    with open("evaled.isophase_SNP.txt", "w") as out_f:
        writer_f = DictWriter(out_f, FIELDS, delimiter="\t")
        writer_f.writeheader()

        debug_count = 0
        if dirs is None:
            dirs = glob.glob("by_loci/*size*/")
        for d1 in dirs:
            # if debug_count > 100: break
            debug_count += 1
            mpileup = Path(d1, "ccs.mpileup")
            mapfile = Path(d1, "fake.mapping.txt")
            vcffile = Path(d1, "phased.partial.vcf")
            nosnp = Path(d1, "phased.partial.NO_SNPS_FOUND")
            if not vcffile.exists():
                assert nosnp.exists()
                logger.info(f"Skipping {d1} because no SNPs found.")
            else:
                logger.info(f"Evaluating {d1}.")
                good_positions, cov_at_pos = get_positions_to_recover(
                    mapfile, mpileup, ki11_snps, min_cov=30
                )  # use lower min cov here becuz a few close cases where BQ filtering lowered cov
                name = d1.split("/")[1]
                eval_isophase(
                    vcffile,
                    ki11_snps,
                    good_positions,
                    cov_at_pos,
                    repeat_by_chrom,
                    ki11_shortread_cov,
                    writer_f,
                    name,
                )

    return ki11_snps


@app.command(name="")
def main(
    vcf_filename: str = typer.Argument(...), out_filename: str = typer.Argument(...)
) -> None:
    brangus(vcf_filename, out_filename)


if __name__ == "__main__":
    typer.run(main)
    # from csv import DictReader

    # main_brangus('combined_97_1modi.vcf', 'evaled.isophase.combined_97_1modi.txt')
    # dirs = ['by_loci/'+r['locus'] for r in DictReader(open('evaled_isophase.demux_hap_count.txt'),delimiter='\t')]
    # main_maize(None, dirs)
