#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
After running collapse, combine with the cluster FL/nFL assignment to get abundance information.
Compliant with SMRTLink/SMRTAnalysis 3.2+ formatting.
May break on older SMRTPortal/SMRTAnalysis 2.x formatting.

cluster_report.csv format in IsoSeq1:

cluster_id,read_id,read_type (FL or NonFL)
i0_ICE_sample03146a|c4,m54026_161015_230744/74318463/31_970_CCS,FL

NOTE: in cluster_report, the sample prefix is i<bin>_ICE_<sample>, but in HQ isoforms is i<bin>_HQ_<sample>!
Must address this :P

cluster_report.csv format in IsoSeq2:

cluster_id,read_id,read_type
cb2729_c2,m54056_171130_193019/8388911/26_1296_CCS,FL


cluster_report.csv format in IsoSeq3:

cluster_id,read_id,read_type
transcript/0,m54056_171130_193019/8388911/ccs,FL


where .group.txt is:

PBfusion.1      HQ_sample0ZPg9hS7|cb7607_c93041/f2p0/1713
PBfusion.2      HQ_sample0ZPg9hS7|cb7607_c16635/f3p0/810,HQ_sample0ZPg9hS7|cb7607_c32934/f2p1/1066

"""


import re
import sys
from collections import Counter
from csv import DictReader, DictWriter
from pathlib import Path
from typing import List, Optional

import typer

from cupcake.__about__ import __version__
from cupcake.logging import cupcake_logger as logger

app = typer.Typer(
    name="get_abundance_post_collapse",
    add_completion=False,
    help="Get abundance/read stat information after running collapse script. Works for Iso-Seq1, 2, and 3 output.",
)


def version_callback(value: bool):
    if value:
        typer.echo(f"cupcake/get_abundance_post_collapse: {__version__}")
        raise typer.Exit()


def get_roi_len(seqid: str):
    # before isoseq3: <movie>/<zmw>/<start>_<end>_CCS
    # for isoseq3: <movie>/<zmw>/ccs
    if seqid.endswith("/ccs"):
        logger.info(
            "WARNING: isoseq3 format detected. Output `length` column will be `NA`."
        )
        return "NA"
    elif not seqid.endswith("_CCS"):
        logger.error(
            "Sequence ID format must be <movie>/<zmw>/<start>_<end>_CCS or <movie>/<zmw>/ccs! Abort!"
        )
        sys.exit(-1)
    s, e, junk = seqid.split("/")[2].split("_")
    return abs(int(s) - int(e))


cluster_rex_sa3 = re.compile(r"(i\d+_ICE_\S+\|c\d+)")
cluster_rex_sa2 = re.compile(r"(i\d+[HL]Q_\S+|c\d+)")
# ex: cb100_c7,m54006_170206_215027/70189550/2661_56_CCS,FL
cluster_rex_tofu2 = re.compile(r"(cb\d+_c\d+)")
cluster_rex_isoseq3 = re.compile(
    r"(transcript/\d+)"
)  # for isoseq3: transcript/0, transcript/1, etc
cluster_rex_isoseq3_mapped = re.compile(
    r"\S+(transcript/\d+)"
)  # for isoseq post SL mapping


def read_group_filename(group_filename: Path, is_cid=True):
    """
    Make the connection between partitioned results and final (ex: PB.1.1)
    The partitioned results could either be ICE cluster (ex: i1_c123) or RoIs

    Group filename format: (SMRTLink/SA3.x, IsoSeq1)
    PB.1.1  i0_HQ_sample03146a|c18072/f2p30/279,i0_HQ_sample03146a|c38435/f2p0/180
    PB.2.1  i2_HQ_sample03146a|c46363/f22p0/2776
    PB.2.2  i2_HQ_sample03146a|c46572/f37p0/2761

    For IsoSeq2:
    PBfusion.1      HQ_sample0ZPg9hS7|cb7607_c93041/f2p0/1713
    PBfusion.2      HQ_sample0ZPg9hS7|cb7607_c16635/f3p0/810,HQ_sample0ZPg9hS7|cb7607_c32934/f2p1/1066

    For IsoSeq3, two possible flavors depending on using conda version or SL Mapping version:
    PB.1.1  transcript/0,transcript/1
    (or)
    PB.1.1  Sample1_RC0_50pM_transcript/0,Sample1_RC0_50pM_transcript/1   <-- in this case we want to strip the prefix out, just keep 'transcript/0'

    Note: LQ/HQ isoforms use the _HQ_ or _LQ_ prefix, but in cluster_report.csv it will be _ICE_
    So we just replace _HQ_ or _LQ_ --> both to be _ICE_

    Returns
    =======
    dict of seq_or_ice_cluster --> collapsed cluster ID (ex: in IsoSeq3 it is 'transcript/0' --> 'PB.1.1')"
    """
    cid_info = {}  # ex: i1 --> c123 --> PB.1.1, or c123 --> PB.1.1

    for line in group_filename.open():
        pbid, members = line.strip().split("\t")
        for cid in members.split(","):
            m = cluster_rex_isoseq3_mapped.match(cid)
            if m is not None:  # IsoSeq3 after mapping, strip it down to 'transcript/0'
                cid = m.group(1)
            else:
                m = cluster_rex_isoseq3.match(cid)
                if m is not None:
                    pass  # pass, nothing to do since the cid is correctly 'transcript/0'
                else:
                    if is_cid:
                        cid = cid.split("/")[0]
                    m = cluster_rex_tofu2.match(cid)
                    if m is not None:
                        cid = m.group(1)
                    else:
                        m = cluster_rex_sa2.match(cid)
                        if m is not None:
                            cid = m.group(1)
                            # replace _HQ_ or _LQ_ with _ICE_ to be compatible with cluster_report.csv later
                            cid = cid.replace("_HQ_", "_ICE_")
                            cid = cid.replace("_LQ_", "_ICE_")
                        elif cluster_rex_sa3.match(cid) is not None:
                            pass  # nothing to do, ID is good
                        elif cid.startswith("HQ_") or cid.startswith("LQ_"):  # isoseq2
                            cid = cid.split("|")[1].split("/")[
                                0
                            ]  # cid = cb7607_c93041, for example
                        else:
                            raise Exception(
                                f"Unrecognized id format {cid} in {group_filename}!"
                            )
            cid_info[cid] = pbid

    return cid_info


def output_read_count_IsoSeq_csv(
    cid_info: List[str],
    csv_filename: str,
    output_filename: str,
    output_mode: Optional[str] = "w",
):
    """
    Turn into read_stats.txt format:
    id \t length \t is_fl \t stat \t pbid
    """
    mapped = {}  # nFL seq -> list of (sample_prefix, cluster) it belongs to

    if output_mode == "w":
        f = open(output_filename, "w")
        f.write("id\tlength\tis_fl\tstat\tpbid\n")
    elif output_mode == "a":
        f = open(output_filename, "a")
    else:
        raise Exception(f"Output mode {output_mode} not valid!")

    unmapped_holder = set()
    for r in DictReader(open(csv_filename), delimiter=","):
        cid = str(r["cluster_id"])
        m = cluster_rex_sa3.match(cid)
        if m is None:
            m = cluster_rex_sa2.match(cid)
            if m is None:
                m = cluster_rex_tofu2.match(cid)
                if m is None:
                    m = cluster_rex_isoseq3.match(cid)
                    if m is not None:
                        cid = m.group(
                            1
                        )  # make sure cid is transcript/0, transcript/1, without any prefix before 'transcript'
                    else:
                        raise Exception(f"cluster_id {cid} is not a valid cluster ID!")

        x = r["read_id"]
        if cid in cid_info:
            # if is FL, must be unique
            if r["read_type"] == "FL":
                pbid, stat = cid_info[cid], "unique"
                f.write(f"{x}\t{get_roi_len(x)}\tY\t{stat}\t{pbid}\n")
            else:  # nonFL could be multi-mapped, must wait and see
                assert r["read_type"] == "NonFL"
                # is only potentially unmapped, add all (movie-restricted) members to unmapped holder
                pbid = cid_info[cid]
                if x not in mapped:
                    mapped[x] = set()
                mapped[x].add(pbid)
        else:
            # unmapped, could be either FL or nFL
            # if FL, can immediately write out since it only appears once in cluster_report.csv
            if r["read_type"] == "FL":
                pbid, stat = "NA", "unmapped"
                f.write(f"{x}\t{get_roi_len(x)}\tY\t{stat}\t{pbid}\n")
            else:  # nonFL could be multi-mapped, so not sure it is truly unmapped, put in holder
                unmapped_holder.add(x)

    # now we can go through the list of mapped/unmapped to see which are uniquely mapped which are not
    for seqid, pbids in mapped.items():
        if len(pbids) == 1:  # unique
            stat = "unique"
        else:
            stat = "ambiguous"
        for pbid in pbids:
            f.write(f"{seqid}\t{get_roi_len(seqid)}\tN\t{stat}\t{pbid}\n")

    unmapped_holder = unmapped_holder.difference(mapped)
    pbid, stat = "NA", "unmapped"
    for x in unmapped_holder:
        f.write(f"{x}\t{get_roi_len(x)}\tN\t{stat}\t{pbid}\n")

    f.close()


def make_abundance_file(
    read_count_filename: str,
    output_filename: str,
    given_total: int = None,
    restricted_movies: List[str] = None,
    write_header_comments: bool = True,
):
    """
    If given_total is not None, use it instead of the total count based on <read_count_filename>
    """
    tally = Counter()  # pbid --> FL count

    reader = DictReader(open(read_count_filename), delimiter="\t")
    for r in reader:
        movie = r["id"].split("/")[0]
        if restricted_movies is None or movie in restricted_movies:
            if r["is_fl"] != "Y":
                logger.info(
                    f'sequence {r["id"]} has `is_fl` field set to `{r["is_fl"]}`. Ignoring this read.'
                )
                continue
            if r["pbid"] == "NA":
                logger.info(
                    f'sequence {r["id"]} has `pbid` field set to `{r["pbid"]}`. Ignoring this read.'
                )
                continue
            assert r["stat"] == "unique"
            tally[r["pbid"]] += 1

    if given_total is not None:
        use_total_fl = given_total
    else:
        use_total_fl = sum(tally.values())

    COUNT_FIELDS = ["pbid", "count_fl", "norm_fl"]
    f = open(output_filename, "w")
    if write_header_comments:
        f.write("#\n")
        f.write("# -----------------\n")
        f.write("# Field explanation\n")
        f.write("# -----------------\n")
        f.write("# count_fl: Number of associated FL reads\n")
        f.write("# norm_fl: count_fl / total number of FL reads, mapped or unmapped\n")
        f.write(f"# Total Number of FL reads: {use_total_fl}\n")
        f.write("#\n")

    writer = DictWriter(f, COUNT_FIELDS, delimiter="\t")
    writer.writeheader()

    keys = list(tally.keys())
    keys.sort(key=lambda x: list(map(int, x.split(".")[1:])))  # sort by PB.1, PB.2....
    for k in keys:
        count_fl = tally[k]
        norm_fl = count_fl * 1.0 / use_total_fl
        rec = {"pbid": k, "count_fl": count_fl, "norm_fl": f"{norm_fl:.4e}"}
        writer.writerow(rec)
    f.close()


def get_abundance_post_collapse(
    group_file: Path,
    cluster_report_csv: Path,
    output_prefix: str,
    restricted_movies: Optional[List[str]] = None,
):
    """

    :param collapse_prefix: collapse prefix filename (must have .group.txt present)
    :param prefix_dict:
    :param output_prefix:
    :param restricted_movies:
    :return:
    """

    if not group_file.exists():
        logger.error(f"File {group_file.name} does not exist. Abort!")
        sys.exit(-1)

    if not cluster_report_csv.exists():
        logger.error(f"File {cluster_report_csv.name} does not exist. Abort!")
        sys.exit(-1)

    cid_info = read_group_filename(group_file, is_cid=True)

    output_read_count_IsoSeq_csv(
        cid_info, cluster_report_csv, f"{output_prefix}.read_stat.txt"
    )
    logger.info(f"Read stat file written to {output_prefix}.read_stat.txt")
    make_abundance_file(
        f"{output_prefix}.read_stat.txt",
        f"{output_prefix}.abundance.txt",
        restricted_movies=restricted_movies,
    )
    logger.info(f"Abundance file written to {output_prefix}.abundance.txt")


@app.command(name="")
def main(
    group_file: Path = typer.Argument(
        ..., help="Group file from collapse_isoforms_by_sam()"
    ),
    cluster_report: Path = typer.Argument(..., help="Cluster CSV report"),
    output_prefix: Optional[str] = typer.Option(
        None,
        help="Name to use for output files.  By default, will use the prefix from the group file",
    ),
    version: Optional[bool] = typer.Option(
        None, "--version", callback=version_callback
    ),
) -> None:
    """Get abundance/read stat information after running collapse script.
    Works for Iso-Seq1, 2, and 3 output."
    """
    if not output_prefix:
        output_prefix = group_file.stem

    get_abundance_post_collapse(group_file, cluster_report, output_prefix)


if __name__ == "__main__":
    typer.run(main)
