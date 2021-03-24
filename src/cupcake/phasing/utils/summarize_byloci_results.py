import glob
from csv import DictReader, DictWriter
from pathlib import Path

import vcf
from Bio import SeqIO
from cupcake.logging import cupcake_logger as logger

FIELDS = ["locus", "size", "num_snp", "num_hap_nopartial", "num_hap_withpartial"]

dirs = [Path(_) for _ in glob.glob("by_loci/*size*")]

with open("summarized.isophase_results.txt", "w") as f:
    writer = DictWriter(f, FIELDS, delimiter="\t")
    writer.writeheader()

    for d in dirs:
        size = 0
        for r in SeqIO.parse(d.joinpath("ccs.fasta").open(), "fasta"):
            size += 1

        rec = {"locus": d, "size": size}

        if d.joinpath(d, "phased.nopartial.NO_SNPS_FOUND").exists():
            rec["num_snp"] = 0
            rec["num_hap_nopartial"] = 0
            rec["num_hap_withpartial"] = 0
        else:
            rec["num_snp"] = len(
                [x for x in vcf.VCFReader(d.joinpath("phased.partial.vcf").open())]
            )
            if d.joinpath("phased.nopartial.NO_HAPS_FOUND").exists():
                rec["num_hap_nopartial"] = 0
                rec["num_hap_withpartial"] = 0
            else:
                file1 = d.joinpath("phased.nopartial.cleaned.human_readable.txt")
                file2 = d.joinpatj("phased.partial.cleaned.human_readable.txt")
                with file1.open() as h1, file2.open() as h2:
                    h1.readline()  # skip header
                    h2.readline()  # skip header
                    rec["num_hap_nopartial"] = len(
                        [r for r in DictReader(h1, delimiter="\t")]
                    )
                    rec["num_hap_withpartial"] = len(
                        [r for r in DictReader(h2, delimiter="\t")]
                    )
        writer.writerow(rec)

logger.info(f"Summarized results of by_loci/<dirs> to {f.name}.")
