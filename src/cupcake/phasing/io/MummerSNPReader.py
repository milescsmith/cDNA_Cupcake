__author__ = "etseng@pacb.com"

"""
For parsing the .snps_files results from running Mummer dnadiff into VCF
https://github.com/mummer4/mummer/blob/master/MANUAL.md#dnadiff

.snps_files file format:
(0) 1-based ref position
(1) base in ref
(2) base in query
(3) 1-based query position
(4) dist to nearest SNP
(5) ???
(6) length of ref
(7) length of query
(8) reading frame of ref (1 is forward, -1 reverse)
(9) reading frame of query
(10) name of ref
(11) name of query

IMPORTANT: Assumes ONLY snps_files and simple indels are in .snps_files. NO LARGE SVS!! This script does not handle them.

There are several diff between .snps_files format and VCF format (see Section 5 for examples)
https://samtools.github.io/hts-specs/VCFv4.2.pdf

for insertions, .snps_files shows ref pos but not the base. VCF expects the ref pos.

46060901        .       C       1679    0       1679    56369465        49753   1       1       000003F 000003F_064
46060901        .       C       1680    0       1680    56369465        49753   1       1       000003F 000003F_064

suppose ref@46060901 is "T"
this in VCF means: pos 46060901, ref: "T", alt: "TCC"

for deletions, .snps_files shows the ref pos but not the query base before it. VCF expects the ref pos before it.

46075044        A       .       15826   1       15826   56369465        49753   1       1       000003F 000003F_064
46075045        G       .       15826   1       15826   56369465        49753   1       1       000003F 000003F_064

suppose ref@46075043 is "T"
this in VCF means: pos 46075043, ref: "TAG", alt: "T"
"""

__VCF_EXAMPLE__ = """
##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT
20      1       .       G       A,T     .       PASS    AF=0.5       GT
"""

from pathlib import Path

import typer
import vcf
from cupcake.logging import cupcake_logger as logger
from cupcake.sequence.SeqReaders import LazyFastaReader
from typer.params import Argument


class SNPRecord(object):
    def __init__(
        self,
        ref_pos,
        query_pos,
        ref_base,
        query_base,
        ref_len,
        query_len,
        ref_frame,
        query_frame,
        ref_name,
        query_name,
    ):
        """
        bases should be 0-based! when printed, will be 1-based.
        """
        self.ref_pos = ref_pos
        self.query_pos = query_pos
        self.ref_base = ref_base
        self.query_base = query_base
        self.ref_len = ref_len
        self.query_len = query_len
        self.ref_strand = "+" if int(ref_frame) == 1 else "-"
        self.query_strand = "+" if int(query_frame) == 1 else "-"
        self.ref_name = ref_name
        self.query_name = query_name

    def __str__(self):
        return f"""
        ref_pos: {self.ref_name}:{self.ref_pos + 1}
        query_pos: {self.query_name}:{self.query_pos + 1}
        ref_base: {self.ref_base}
        query_base: {self.query_base}
        """


class SNPReader(object):
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename)

    def __iter__(self):
        return self

    def __next__(self):
        cur = self.f.tell()
        line = self.f.readline()
        if self.f.tell() == cur:
            raise StopIteration
        return self.parseLine(line)

    def parseLine(self, line):
        raw = line.strip().split("\t")
        if len(raw) != 12:
            raise Exception(
                f"Expected to have 12 cols in MUMMER SNP record \
            but saw only {len(raw)}, abort! Line was: {line}"
            )
        return SNPRecord(
            ref_pos=int(raw[0]) - 1,
            query_pos=int(raw[3]) - 1,
            ref_base=raw[1],
            query_base=raw[2],
            ref_len=int(raw[6]),
            query_len=int(raw[7]),
            ref_frame=int(raw[8]),
            query_frame=int(raw[9]),
            ref_name=raw[10],
            query_name=raw[11],
        )


app = typer.Typer(
    name="cupcake.phasing.io.MummerSNPReader",
    help="Process one or more .snps_files files from dnadiff to VCF format.",
)


def write_snp_to_vcf(
    snp_filename: Path,
    vcf_filename: Path,
    genome_filename: Path,
    genome_d: LazyFastaReader = None,
) -> None:
    # read the genome is genome_d is not given
    if genome_d is None:
        genome_d = LazyFastaReader(genome_filename)

    # read the first SNP record so we know the query name
    snp_reader = SNPReader(snp_filename)
    snp_rec = next(snp_reader)
    sample_name = snp_rec.query_name
    cur_recs = [snp_rec]
    genome_rec = genome_d[snp_rec.ref_name]

    with open("template.vcf", "w+") as f:
        f.write(f"{__VCF_EXAMPLE__}\n")
        reader = vcf.VCFReader(f)
        reader.samples = [sample_name]
        f_vcf = vcf.Writer(vcf_filename.open("w"), reader)

        for r1 in snp_reader:
            if r1.ref_pos == cur_recs[-1].ref_pos:  # multi-nt insertion, keep recording
                cur_recs.append(r1)
            elif (
                r1.query_base == "." and cur_recs[-1].query_base == "."
            ):  # multi-nt deletion, keep recording
                cur_recs.append(r1)
            else:  # time to write out the current set of records
                # multiple records mean it could be:
                # 1. multi-nucleotide insertions
                # 2. multi-nucleotide deletions

                if (
                    len(cur_recs) == 1
                    and cur_recs[0].ref_base != "."
                    and cur_recs[0].query_base != "."
                ):  # just a SNP record
                    pos = cur_recs[0].ref_pos
                    ref_base = cur_recs[0].ref_base
                    alt_base = cur_recs[0].query_base
                elif cur_recs[0].ref_base == ".":
                    # is a single or multi-nt insertions, must retrieve ref base from genome
                    # ex: in out.snps_files it is . --> ATG
                    # in VCF it should be T --> TATG (meaning insertion of ATG)
                    pos = cur_recs[0].ref_pos
                    ref_base = genome_rec[cur_recs[0].ref_pos]
                    alt_base = ref_base + "".join(r.query_base for r in cur_recs)
                else:
                    # is a single multi-nt deletions, we need to get one more ref base before the first deletion
                    # ex: in out.snps_files it is GGG --> deletion
                    # in VCF it should be TGGG --> T (meaning deletion of GGG)
                    pos = cur_recs[0].ref_pos - 1
                    ref_base_prev = genome_rec[pos]
                    ref_base = ref_base_prev + "".join(r.ref_base for r in cur_recs)
                    alt_base = ref_base_prev

                rec = vcf.model._Record(
                    CHROM=snp_rec.ref_name,
                    POS=pos + 1,
                    ID=".",
                    REF=ref_base,
                    ALT=[vcf.model._Substitution(alt_base)],
                    QUAL=".",
                    FILTER="PASS",
                    INFO={"AF": 0.5},
                    FORMAT="GT",
                    sample_indexes=None,
                )
                samp_ft = vcf.model.make_calldata_tuple(["GT"])
                rec.samples.append(vcf.model._Call(rec, sample_name, samp_ft(*["0|1"])))
                f_vcf.write_record(rec)
                if r1.ref_name != cur_recs[0].ref_name:
                    genome_rec = genome_d[r1.ref_name]
                cur_recs = [r1]


@app.command(name="")
def main(
    snps_filename: str = typer.Argument(
        ..., help="Filename containing list of .snps_files to process."
    ),
    genome_filename: str = typer.Argument(
        ...,
        help="Genome fasta. Chromosome IDs must agree with .snps_files files!",
    ),
):
    snps_filename = Path(snps_filename)
    genome_filename = Path(genome_filename)

    snps_files = []
    # sanity checking of input files
    for filename in snps_filename.open():
        if not filename.suffix(".snps"):
            raise FileNotFoundError(
                f"Input files listed in {snps_filename} must end with .snps_files!"
            )
        if not filename.exists():
            raise FileNotFoundError(f"{filename} does not exist! Abort.")
        snps_files.append(filename)

    if not genome_filename.exists():
        raise FileNotFoundError(f"Genome file {genome_filename} does not exist!")

    logger.info(f"Reading genome file {genome_filename}....")
    genome_d = LazyFastaReader(genome_filename)

    # quick checking if the genome chromosomes have the |arrow|arrow style suffix, if they do, process it
    keys = list(genome_d.keys())
    for k in keys:
        k2 = k.split("|")[0]
        if k2 != k and k2 not in keys:
            genome_d.d[k2] = genome_d.d[k]
            logger.info(
                f"Detected | string in chromosome ID, stripping {k} to {k2}...."
            )
    logger.info("Finished reading genome.")

    for snp_file in snps_files:
        assert snp_file.suffix(".snps")
        vcf_file = snp_file.with_suffix(".vcf")
        logger.info(f"Processing {snp_file} --> {vcf_file}")
        write_snp_to_vcf(snp_file, vcf_file, genome_filename, genome_d)


if __name__ == "__main__":
    typer.run(main)()
