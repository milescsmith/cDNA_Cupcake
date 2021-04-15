import glob
from collections import Counter, defaultdict
from enum import Enum
from pathlib import Path
from typing import Optional

import typer
import vcfpy

from cupcake.logging import cupcake_logger as logger


class vcf_phasing(str, Enum):
    phased_partial_vcf = "phased.partial.vcf"
    phased_partial_cleaned_vcf = "phased.partial.cleaned.vcf"
    phased_nopartial_vcf = "phased.nopartial.vcf"
    phased_nopartial_cleaned_vcf = "phased.nopartial.cleaned.vcf"


app = typer.Typer(name="cupcake.phasing.utils.collect_all_vcf")


def collect_all_vcf(
    dirs: str,
    vcf_filename: str = "phased.partial.vcf",
    output: str = "IsoSeq_IsoPhase.vcf",
) -> None:
    no_snp_found_filename = Path(f"{Path(vcf_filename).stem}.NO_SNPS_FOUND")
    snps_by_chrom = defaultdict(lambda: [])

    reader = None

    for d in dirs:
        filename = Path(d, vcf_filename)
        if not filename.exists():
            if not no_snp_found_filename.exists():
                logger.info("VCF file {filename} does not exist. Skipping.")
            continue
        with open(filename) as rf:
            reader = vcfpy.Reader(rf)

            for r in reader:
                c = Counter()  # genotype -> count
                for x in r.samples:
                    if x.data.GT.count("|") == 0:
                        c[x.data.GT] += x.data.HQ
                    else:
                        for i, gt in enumerate(x.data.GT.split("|")):
                            c[gt] += x.data.HQ[i]
                c_keys = c.keys()
                genotype = "|".join(str(k) for k in c_keys)
                counts = ",".join(str(c[k]) for k in c_keys)
                r.samples = [
                    vcfpy.Call(
                        r,
                        "SAMPLE",
                        vcfpy.OrderedDict([("GT", genotype), ("HQ", counts)]),
                    )
                ]
                snps_by_chrom[r.CHROM].append((r.POS, r))

    keys = list(snps_by_chrom.keys())
    keys.sort()

    if reader is not None:
        reader.samples = ["SAMPLE"]
        with open(output, "w") as f:
            f = vcfpy.Writer(f, reader)
            for k in keys:
                v = snps_by_chrom[k]
                v.sort(key=lambda x: x[0])
                for _, rec in v:
                    f.write_record(rec)
        print("Output written to:", output)


@app.command(name="")
def main(
    directory: str = typer.Option(
        "by_loci/",
        "--dir",
        "-d",
        help="Directory containing the subdirs of IsoPhase",
    ),
    output: str = typer.Option(
        "IsoSeq_IsoPhase.vcf",
        "--output",
        "-o",
        help="Output VCF filename (default: IsoSeq_IsoPhase.vcf)",
    ),
    vcf: vcf_phasing = typer.Option(
        vcf_phasing.phased_partial_cleaned_vcf,
        help="VCF to use per directory (default: phased.partial.cleaned.vcf)",
    ),
    select_dirs: Optional[str] = typer.Option(
        None,
        "--select_dirs",
        "-s",
        help="Comma separate list of directories to tally - if this is used, <dir> is ignored",
    ),
) -> None:

    if select_dirs is None:
        dirs = glob.glob(f"{directory}/*size*")
    else:
        dirs = select_dirs.split(",")
        for d in dirs:
            if not Path(d).is_dir() or not Path(d).exists():
                raise FileNotFoundError(
                    f"{d} is not a directory or does not exist! Abort!"
                )

    collect_all_vcf(dirs, vcf, output)


if __name__ == "__main__":
    typer.run(main)
