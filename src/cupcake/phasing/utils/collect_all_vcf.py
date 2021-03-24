import glob
import sys
from collections import Counter, defaultdict
from pathlib import Path

import vcf
from cupcake.logging import cupcake_logger as logger


def collect_all_vcf(
    dirs: str,
    vcf_filename: str = "phased.partial.vcf",
    output: str = "IsoSeq_IsoPhase.vcf",
) -> None:
    no_snp_found_filename = Path(f"{Path(vcf_filename).stem}.NO_SNPS_FOUND")
    snps_by_chrom = defaultdict(lambda: [])

    samp_ft = vcf.model.make_calldata_tuple(["GT", "HQ"])
    reader = None

    for d in dirs:
        filename = Path(d, vcf_filename)
        if not filename.exists():
            if not no_snp_found_filename.exists():
                logger.info("VCF file {filename} does not exist. Skipping.")
            continue
        with open(filename) as rf:
            reader = vcf.VCFReader(rf)

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
                r.samples = [vcf.model._Call(r, "SAMPLE", samp_ft(*[genotype, counts]))]
                snps_by_chrom[r.CHROM].append((r.POS, r))

    keys = list(snps_by_chrom.keys())
    keys.sort()

    if reader is not None:
        reader.samples = ["SAMPLE"]
        with open(output, "w") as f:
            f = vcf.Writer(f, reader)
            for k in keys:
                v = snps_by_chrom[k]
                v.sort(key=lambda x: x[0])
                for pos, rec in v:
                    f.write_record(rec)
        print("Output written to:", output)


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument(
        "-d",
        "--dir",
        default="by_loci/",
        help="Directory containing the subdirs of IsoPhase (default: by_loci/)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="IsoSeq_IsoPhase.vcf",
        help="Output VCF filename (default: IsoSeq_IsoPhase.vcf)",
    )
    parser.add_argument(
        "--vcf",
        default="phased.partial.cleaned.vcf",
        choices=[
            "phased.partial.vcf",
            "phased.partial.cleaned.vcf",
            "phased.nopartial.vcf",
            "phased.nopartial.cleaned.vcf",
        ],
        help="VCF to use per directory (default: phased.partial.cleaned.vcf)",
    )
    parser.add_argument(
        "-s",
        "--select_dirs",
        default=None,
        help="Comma separate list of directories to tally - if this is used, <dir> is ignored",
    )

    args = parser.parse_args()

    if args.select_dirs is None:
        dirs = glob.glob(f"{args.dir}/*size*")
    else:
        dirs = args.select_dirs.split(",")
        for d in dirs:
            if not Path(d).is_dir() or not Path(d).exists():
                logger.error(f"{d} is not a directory or does not exist! Abort!")
                sys.exit(-1)

    collect_all_vcf(dirs, args.vcf, args.output)


if __name__ == "__main__":
    main()
