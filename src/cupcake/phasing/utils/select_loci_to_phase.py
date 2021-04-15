__author__ = "etseng@pacb.com"

"""
For Iso-Phase, selecting loci that has sufficient FL coverage to phase.

INPUT:
   -- FLNC (fastq)
   -- Unique isoforms and genes (GFF), using PB.X.Y to denote gene vs isoforms relationship
   -- FLNC association to isoforms, (.read_stat.txt)
   -- Reference genome (fasta)
"""

import re
import sys
from collections import defaultdict, namedtuple
from csv import DictReader
from pathlib import Path

import typer
from Bio import SeqIO
from bx.intervals.cluster import ClusterTree

from cupcake.logging import cupcake_logger as logger
from cupcake.sequence.GFF import collapseGFFReader
from cupcake.sequence.SeqReaders import LazyFastqReader

rex_flnc = re.compile(r"(m\S+_\d+_\d+\/\d+)\/ccs")  # Sequel/Iso-Seq 3 format
rex_flnc2 = re.compile(
    r"(m\d+_\d+_\d+_\w+_s1_p0\/\d+)\/ccs"
)  # ex: m160920_210440_42165_c101101052550000001823258304261787_s1_p0/83826/ccs
rex_flnc3 = re.compile(
    r"(m\d+_\d+_\d+_\w+_s1_p0\/\d+)\/\d+_\d+_CCS"
)  # ex: m160918_184656_42165_c101101052550000001823258304261782_s1_p0/121477/5106_57_CCS
rex_pbid = re.compile(r"(PB.\d+).\d+")


extra_bp_around_junctions = 50  # get this much around junctions to be safe AND to not screw up GMAP who doesn't like microintrons....
__padding_before_after__ = 10  # get this much before and after the start


app = typer.Typer(name="cupcake.phasing.utils.select_loci_to_phase")


def read_flnc_fastq(flnc_filename):
    """
    Read FLNC fastq into a dict of zmw --> lazy file pointer
    """
    flnc_fastq_d = LazyFastqReader(flnc_filename)
    rich_zmws = set()

    for k in list(flnc_fastq_d.keys()):
        m = rex_flnc.match(k)
        if m is None:
            m = rex_flnc2.match(k)
            if m is None:
                m = rex_flnc3.match(k)
                if m is None:
                    raise Exception(
                        f"Expected FLNC id format is <movie>/<zmw>/ccs! Instead saw: {k}"
                    )
        zmw = m.group(1)
        flnc_fastq_d.d[zmw] = flnc_fastq_d.d[k]
        rich_zmws.add(zmw)

    return flnc_fastq_d, rich_zmws


def read_read_stat(stat_filename, rich_zmws):
    """
    Read .read_stat.txt file
    :return: tally_by_loci -- dict  of {locus -- list of (isoform, zmw)}
             poor_zmws_not_in_rich --- list of ZMWs that were in .read_stat but missing in FLNC FASTQ. ONly happens if used diff CCS to get FASTQ.
    """
    required_fields = ["pbid", "id", "is_fl", "stat"]
    tally_by_loci = defaultdict(
        lambda: []
    )  # PB.1 --> [(PB.1.1, zmw1), (PB.1.1, zmw2), (PB.1.2, zmw3)...]
    poor_zmws_not_in_rich = set()
    reader = DictReader(open(stat_filename), delimiter="\t")
    if any(x not in reader.fieldnames for x in required_fields):
        raise Exception(
            f"Expected fields `pbid`, `is_fl`, and `stat` in {stat_filename} but only saw: {reader.fieldnames}"
        )

    for r in reader:
        if r["is_fl"] == "Y" and r["stat"] == "unique":
            m = rex_pbid.match(r["pbid"])
            if m is None:
                raise Exception(f"Expected PBID format PB.X.Y but saw {r['pbid']}")
            locus = m.group(1)  # ex: PB.1
            m = rex_flnc.match(r["id"])
            if m is None:
                m = rex_flnc2.match(r["id"])
                if m is None:
                    m = rex_flnc3.match(r["id"])
                    if m is None:
                        raise Exception(
                            f"Expected FLNC id format is <movie>/<zmw>/ccs! Instead saw: {r['id']}"
                        )
            zmw = m.group(1)
            if zmw in rich_zmws:
                tally_by_loci[locus].append((r["pbid"], zmw))
            else:
                poor_zmws_not_in_rich.add(zmw)

    return tally_by_loci, poor_zmws_not_in_rich


LocusInfo = namedtuple("LocusInfo", ["chrom", "strand", "regions", "isoforms"])


def read_GFF(gff_filename):
    """
    Read a GFF filename and get the gene regions

    :return: dict of (PB.X) --> LocusInfo
    """
    gff_info = {}  # loci --> LocusInfo
    tmp = {}  # loci PB.X --> list of GFF records for PB.X.Y

    for r in collapseGFFReader(gff_filename):
        m = rex_pbid.match(r.seqid)
        if m is None:
            raise Exception(f"Expected PBID format PB.X.Y but saw {r.seqid}")
        locus = m.group(1)  # ex: PB.1
        if locus not in tmp:
            tmp[locus] = [r]
            gff_info[locus] = LocusInfo(
                chrom=r.chr, strand=r.strand, regions=None, isoforms=None
            )
        else:
            if gff_info[locus].chrom != r.chr:
                logger.warning(
                    f"WARNING: Expected {r.seqid} to be on {gff_info[locus].chrom} but saw {r.chr}. Could be minimap2 multi-mapping inconsistency for repetitive genes. Check later.\n"
                )
            tmp[locus].append(r)

    # now figure out the exonic regions for each gene PB.X
    for locus, records in tmp.items():
        c = ClusterTree(0, 0)
        for r in records:
            for e in r.ref_exons:
                c.insert(
                    max(0, e.start - extra_bp_around_junctions),
                    e.end + extra_bp_around_junctions,
                    1,
                )

        regions = [(a, b) for (a, b, junk) in c.getregions()]
        regions[0] = (max(0, regions[0][0] - __padding_before_after__), regions[0][1])
        regions[-1] = (
            max(0, regions[-1][0]),
            regions[-1][1] + __padding_before_after__,
        )
        gff_info[locus] = LocusInfo(
            chrom=gff_info[locus].chrom,
            strand=gff_info[locus].strand,
            regions=regions,
            isoforms=[r.seqid for r in records],
        )

    return gff_info


def make_fake_genome(genome_d, gff_info, locus, output_prefix, output_name):

    chrom = gff_info[locus].chrom
    regions = gff_info[locus].regions

    with open(f"{output_prefix}.fasta", "w") as f:
        f.write(">" + output_name + "\n")
        for s, e in regions:
            f.write(str(genome_d[chrom][s:e].seq))
        f.write("\n")
        f.close()

    # for mapping, write <0-based index on fake genome>, <ref chrom>, <0-based index on ref genome>
    with open(f"{output_prefix}.mapping.txt", "w") as f:
        i = 0
        for s, e in regions:
            for j in range(s, e):
                f.write(f"{i},{chrom},{j}\n")
                i += 1

    with open(f"{output_prefix}.pbids.txt", "w") as f:
        f.write("\n".join(gff_info[locus].isoforms) + "\n")

    logger.info(
        f"Output written to {output_prefix}.fasta, {output_prefix}.mapping.txt, {output_prefix}.pbids.txt.",
    )


def select_loci_to_phase(
    genome_dict,
    gff_filename,
    stat_filename,
    flnc_filename,
    coverage,
):

    logger.info(
        "Reading FLNC file...",
    )
    flnc_fastq_d, rich_zmws = read_flnc_fastq(flnc_filename)

    logger.info(
        "Reading read_stat...",
    )
    tally_by_loci, poor_zmws_not_in_rich = read_read_stat(stat_filename, rich_zmws)

    logger.info(
        "Reading GFF file...",
    )
    gff_info = read_GFF(gff_filename)

    # find all gene loci that has at least X FLNC coverage
    cand_loci = [k for k in tally_by_loci if len(tally_by_loci[k]) >= coverage]
    logger.info(
        f"Total {len(tally_by_loci)} loci read. {len(cand_loci)} has >= {coverage} coverage.",
    )
    for locus in cand_loci:
        if locus not in gff_info:
            logger.warning(
                f"WARNING: {locus} skipped because not in GFF info (probably filtered out later).\n"
            )
            continue

        logger.info(f"making {locus}")
        d2 = Path(f"by_loci/{locus}_size{len(tally_by_loci[locus])}")
        d2.mkdir()

        chrom_len = {k: len(v) for k, v in genome_dict.items()}

        with d2.joinpath("config").open("w") as f:
            ref_start = max(0, gff_info[locus].regions[0][0])  # ref start must be >= 0
            # ref end must be at most chrom length
            ref_end = min(
                chrom_len[gff_info[locus].chrom], gff_info[locus].regions[-1][-1]
            )

            # _chr, _strand, _start, _end = gff_info[locus]
            f.write(f"pbid={locus}\n")
            f.write(f"ref_chr={gff_info[locus].chrom}\n")
            f.write(f"ref_strand={gff_info[locus].strand}\n")
            f.write(f"ref_start={ref_start}\n")
            f.write(f"ref_end={ref_end}\n")

        make_fake_genome(genome_dict, gff_info, locus, f"{d2}/fake", f"fake_{locus}")

        # write ccs.fastq
        with d2.joinpath("ccs.fastq").open("w") as f1, d2.joinpath("ccs.fasta").open(
            "w"
        ) as f2, d2.joinpath("fake.read_stat.txt").open("w") as h:
            h.write("id\tlength\tis_fl\tstat\tpbid\n")
            for pbid, zmw in tally_by_loci[locus]:
                rec = flnc_fastq_d[zmw]
                SeqIO.write(rec, f1, "fastq")
                SeqIO.write(rec, f2, "fasta")
                h.write(f"{zmw}\t{len(rec.seq)}\tY\tunique\t{pbid}\n")


@app.command(name="")
def main(
    genome_fasta: str = typer.Argument(..., help="Reference genome fasta"),
    flnc_filename: str = typer.Argument(..., help="FLNC fastq file"),
    gff_filename: str = typer.Argument(
        ..., help="GFF file of transcripts, IDs must be PB.X.Y"
    ),
    stat_filename: str = typer.Argument(
        ..., help="Tab-delimited read stat file linking FLNC to PB.X.Y"
    ),
    coverage: int = typer.Option(
        40,
        "--coverage",
        "-c",
        help="Minimum FLNC coverage required (default: 40)",
    ),
) -> None:
    if Path("by_loci").exists() and Path("by_loci").is_dir():
        logger.error("Directory by_loci/ already exists. Delete before running!")
        sys.exit(-1)

    if not Path(genome_fasta).exists():
        logger.error(f"Cannot find genome FASTA {genome_fasta}. Abort!")
        sys.exit(-1)

    if not Path(flnc_filename).exists():
        logger.error(f"Cannot find FLNC file {flnc_filename}. Abort!")
        sys.exit(-1)

    if not Path(gff_filename).exists():
        logger.error(f"Cannot find GFF file {gff_filename}. Abort!")
        sys.exit(-1)

    if not Path(stat_filename).exists():
        logger.error(f"Cannot find Stat file {stat_filename}. Abort!")
        sys.exit(-1)

    logger.info(f"Reading genome fasta {genome_fasta}...")
    genome_d = SeqIO.to_dict(SeqIO.parse(open(genome_fasta), "fasta"))

    select_loci_to_phase(genome_d, gff_filename, stat_filename, flnc_filename, coverage)


if __name__ == "__main__":
    typer.run(main)
