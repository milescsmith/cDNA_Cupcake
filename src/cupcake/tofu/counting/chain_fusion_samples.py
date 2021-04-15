#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Similar to chain_samples.py, but specifically for chaining fusions.

Difference with regular chaining:

1. A fusion consists of several transcript records, denoted by PB.X.1, PB.X.2, ... where PB.X is one fusion
   (must always consider possibility of a fusion consisting of 2+ loci)

2. For two fusions from two samples to be "equivalent", every loci (.1, .2, ..) must agree

"""

import shutil
from collections import defaultdict
from csv import DictReader
from enum import Enum
from pathlib import Path
from typing import List, Optional, Union

import typer
from Bio import SeqIO

import cupcake.sequence.GFF as GFF
from cupcake.logging import cupcake_logger as logger
from cupcake.tofu.counting import combine_abundance_across_samples as sp
from cupcake.tofu.counting.chain_samples import read_config, read_count_info

app = typer.Typer(name="cupcake.tofu.counting.chain_fusion_sample")


class fl_fields(str, Enum):
    norm_fl = "norm_fl"
    count_fl = "count_fl"


def sample_sanity_check(
    group_filename: Union[str, Path],
    gff_filename: Union[str, Path],
    count_filename: Union[str, Path],
    fastq_filename: Optional[Union[str, Path]] = None,
) -> None:
    """
    Double check that the formats are expected and all PBIDs are concordant across the files
    :return: raise Exception if sanity check failed
    """

    logger.info(
        f"Sanity checking. Retrieving PBIDs from {group_filename},{gff_filename},{count_filename}..."
    )
    ids1 = [line.strip().split()[0] for line in open(group_filename)]
    ids2 = [fusion_id for fusion_id, rs in GFF.collapseGFFFusionReader(gff_filename)]
    with open(count_filename) as f:
        for _ in range(14):
            f.readline()  # just through the header
        ids3 = [r["pbid"] for r in DictReader(f, delimiter="\t")]
        if len({ids2}.difference(ids1)) > 0 or len({ids2}.difference(ids3)) > 0:
            raise Exception(
                f"Sanity check failed! Please make sure the PBIDs listed in {gff_filename} are also in {gff_filename} and {count_filename}"
            )

    if fastq_filename is not None:
        ids4 = [r.id.split("|")[0] for r in SeqIO.parse(fastq_filename, "fastq")]
        if len({ids2}.difference(ids4)) > 0:
            raise Exception(
                f"Sanity check failed! Please make sure the PBIDs listed in {gff_filename} are also in {fastq_filename}"
            )


def chain_fusion_samples(
    dirs: List[str],
    names: List[str],
    group_filename: Union[str, Path],
    gff_filename: Union[str, Path],
    count_filename: Union[str, Path],
    field_to_use: str = "count_fl",
    fuzzy_junction: int = 0,
    fastq_filename: Optional[Union[str, Path]] = None,
) -> None:

    for d in dirs.values():
        sample_sanity_check(
            Path(d, group_filename),
            Path(d, gff_filename),
            Path(d, count_filename),
            Path(d, fastq_filename) if fastq_filename is not None else None,
        )

    _, count_info = read_count_info(count_filename, dirs, field_to_use)

    # some names may already start with "tmp_" which means they are intermediate results that have already been chained
    # find the first non "tmp_" and start from there
    if names[0].startswith("tmp_"):
        chain = []
        for start_i, name in enumerate(names):
            if name.startswith("tmp_"):
                chain.append(name[4:])
            else:
                break
        # start_i, name now points at the first "non-tmp" sample
        # we want to go to the last tmp_ sample and read it
        name = names[start_i - 1][4:]  # this is the last tmp_ sample, let's read it
        o = sp.MegaPBTreeFusion(
            gff_filename=f"tmp_{name}.gff",
            group_filename=f"tmp_{name}.group.txt",
            self_prefix=f"tmp_{name}",
            internal_fuzzy_max_dist=fuzzy_junction,
            fastq_filename=f"tmp_{name}.rep.fq" if fastq_filename is not None else None,
        )
        # chain.append(name) # no need, already done above
    else:  # everything is new, start fresh
        name = names[0]
        d = dirs[name]
        chain = [name]
        o = sp.MegaPBTreeFusion(
            gff_filename=Path(d, gff_filename),
            group_filename=Path(d, group_filename),
            self_prefix=name,
            internal_fuzzy_max_dist=fuzzy_junction,
            fastq_filename=Path(d, fastq_filename)
            if fastq_filename is not None
            else None,
        )
        start_i = 1

    for name in names[start_i:]:
        assert not name.startswith("tmp_")
        d = dirs[name]
        o.add_sample(
            gff_filename=Path(d, gff_filename),
            group_filename=Path(d, group_filename),
            sample_prefix=name,
            output_prefix=f"tmp_{name}",
            fastq_filename=Path(d, fastq_filename)
            if fastq_filename is not None
            else None,
        )
        o = sp.MegaPBTreeFusion(
            gff_filename=f"tmp_{name}.gff",
            group_filename=f"tmp_{name}.group.txt",
            self_prefix=f"tmp_{name}",
            internal_fuzzy_max_dist=fuzzy_junction,
            fastq_filename=f"tmp_{name}.rep.fq" if fastq_filename is not None else None,
        )
        chain.append(name)

    # now recursively chain back by looking at mega_info.txt!!!
    d = {}  # ex: (tmp_1009, PB.1.1) --> mega info dict
    for c in chain[1:]:
        for r in DictReader(open(f"tmp_{c}.mega_info.txt"), delimiter="\t"):
            d[f"tmp_{c}", r["pbid"]] = r

    with open("all_samples.chained_ids.txt", "w") as f1, open(
        "all_samples.chained_count.txt", "w"
    ) as f2:
        f1.write("superPBID")
        f2.write("superPBID")
        for c in chain:
            f1.write(f"	{c}")
            f2.write(f"	{c}")
        f1.write("\n")
        f2.write("\n")

        reader = DictReader(
            Path(f"tmp_{chain[-1]}.mega_info.txt").open(), delimiter="\t"
        )
        for r in reader:
            saw_NA = False
            r0 = r
            answer = defaultdict(lambda: "NA")  # ex: 1009 --> PB.1.1
            answer2 = defaultdict(lambda: "NA")  # ex: 1009 --> count
            answer[chain[-1]] = r[chain[-1]]
            if r[chain[-1]] != "NA":
                answer2[chain[-1]] = count_info[chain[-1], answer[chain[-1]]]
            for c in chain[::-1][
                1:-1
            ]:  # the first sample does not have tmp_, because it's not a chain
                if r[f"tmp_{c}"] == "NA":
                    saw_NA = True
                    break
                else:
                    r2 = d[f"tmp_{c}", r[f"tmp_{c}"]]
                    answer[c] = r2[c]
                    if answer[c] != "NA":
                        answer2[c] = count_info[c, answer[c]]
                    r = r2
            if not saw_NA:
                answer[chain[0]] = r[chain[0]]
                if answer[chain[0]] != "NA":
                    answer2[chain[0]] = count_info[chain[0], answer[chain[0]]]
            f1.write(r0["pbid"])
            f2.write(r0["pbid"])
            for c in chain:
                f1.write(f"	{answer[c]}")  # each tissue still share the same PB id
                f2.write(f"	{str(answer2[c])}")
            f1.write("\n")
            f2.write("\n")

    shutil.copyfile(f"tmp_{chain[-1]}.gff", "all_samples.chained.gff")
    if fastq_filename is not None:
        shutil.copyfile(f"tmp_{chain[-1]}.rep.fq", "all_samples.chained.rep.fq")

    logger.info("Chained output written to:")
    logger.info("all_samples.chained.gff")
    logger.info(f1.name)
    logger.info(f2.name)
    if fastq_filename is not None:
        logger.info("all_samples.chained.rep.fq")


@app.command(name="")
def main(
    config_file: str = typer.Argument(...),
    field_to_use: fl_fields = typer.Option(
        fl_fields.count_fl,
        show_default=False,
        help="Which count field to use for chained sample",
    ),
    fuzzy_junction: int = typer.Option(
        5,
        show_default=False,
        help="Max allowed distance in junction to be considered identical",
    ),
) -> None:
    (
        sample_dirs,
        sample_names,
        group_filename,
        gff_filename,
        count_filename,
        fastq_filename,
    ) = read_config(config_file)

    chain_fusion_samples(
        sample_dirs,
        sample_names,
        group_filename,
        gff_filename,
        count_filename,
        field_to_use,
        fuzzy_junction,
        fastq_filename,
    )


if __name__ == "__main__":
    typer.run(main)
