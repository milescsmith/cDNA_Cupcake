#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import shutil
import sys
from collections import OrderedDict, defaultdict
from csv import DictReader, DictWriter
from enum import Enum
from multiprocessing import Process
from pathlib import Path
from typing import List, Tuple, Union

import typer
from Bio import SeqIO
from bx.intervals.cluster import ClusterTree

from cupcake.logging import cupcake_logger as logger
from cupcake.sequence import GFF
from cupcake.tofu.counting import combine_abundance_across_samples as sp

app = typer.Typer(name="cupcake.tofu.counting.chain_samples")


class fl_fields(str, Enum):
    norm_fl = "norm_fl"
    count_fl = "count_fl"


def sample_sanity_check(
    group_filename, gff_filename, count_filename, fastq_filename=None
) -> None:
    """
    Double check that the formats are expected and all PBIDs are concordant across the files
    :return: raise Exception if sanity check failed
    """
    logger.info(
        f"Sanity checking. Retrieving PBIDs from {group_filename},{gff_filename},{count_filename}..."
    )
    ids1 = [line.strip().split()[0] for line in open(group_filename)]
    ids2 = [r.seqid for r in GFF.collapseGFFReader(gff_filename)]
    f = open(count_filename)
    while True:
        # advance through the headers which start with #
        cur = f.tell()
        if (
            not f.readline().startswith("#") or f.tell() == cur
        ):  # first non-# seen or EOF
            f.seek(cur)
            break
    ids3 = [r["pbid"] for r in DictReader(f, delimiter="\t")]
    if len(set(ids2).difference(ids1)) > 0 or len(set(ids2).difference(ids3)) > 0:
        raise Exception(
            f"Sanity check failed! Please make sure the PBIDs listed in {gff_filename} are also in {group_filename} and {count_filename}"
        )

    if fastq_filename is not None:
        ids4 = [r.id.split("|")[0] for r in SeqIO.parse(open(fastq_filename), "fastq")]
        if len(set(ids2).difference(ids4)) > 0:
            raise Exception(
                f"Sanity check failed! Please make sure the PBIDs listed in {gff_filename} are also in {fastq_filename}"
            )


def read_config(filename):
    """
    tmpSAMPLE=<name>;<path>
    SAMPLE=<name>;<path>

    must also have
    GROUP_FILENAME=
    GFF_FILENAME=
    COUNT_FILENAME=

    optional:
    FASTQ_FILENAME=
    """
    sample_dirs = {}
    sample_names = []
    group_filename, gff_filename, count_filename = None, None, None
    fastq_filename = None

    no_more_tmp = False

    with open(filename) as f:
        for line in f:
            if line.startswith("tmpSAMPLE="):
                if no_more_tmp:
                    logger.error(
                        "Cannot have tmp_ samples after non-tmp_ samples! Abort!"
                    )
                    sys.exit(-1)
                name, path = line.strip()[len("tmpSAMPLE=") :].split(";")
                if name.startswith("tmp_"):
                    logger.error(
                        f"Sample names are not allowed to start with tmp_! "
                        f"Please change {name} to something else."
                    )
                    sys.exit(-1)
                sample_dirs[name] = Path(path).resolve()
                sample_names.append(f"tmp_{name}")
            elif line.startswith("SAMPLE="):
                no_more_tmp = True
                name, path = line.strip()[len("SAMPLE=") :].split(";")
                if name.startswith("tmp_"):
                    logger.error(
                        f"Sample names are not allowed to start with tmp_! "
                        f"Please change {name} to something else."
                    )
                    sys.exit(-1)
                sample_dirs[name] = Path(path).resolve()
                sample_names.append(name)
            elif line.startswith("GROUP_FILENAME="):
                group_filename = line.strip()[len("GROUP_FILENAME=") :]
            elif line.startswith("GFF_FILENAME="):
                gff_filename = line.strip()[len("GFF_FILENAME=") :]
            elif line.startswith("COUNT_FILENAME="):
                count_filename = line.strip()[len("COUNT_FILENAME=") :]
            elif line.startswith("FASTQ_FILENAME="):
                fastq_filename = line.strip()[len("FASTQ_FILENAME=") :]

    if group_filename is None:
        raise Exception(
            f"Expected GROUP_FILENAME= but not in config file {filename}! Abort."
        )
    if count_filename is None:
        raise Exception(
            f"Expected COUNT_FILENAME= but not in config file {filename}! Abort."
        )
    if gff_filename is None:
        raise Exception(
            f"Expected GFF_FILENAME= but not in config file {filename}! Abort."
        )

    if len(sample_names) == 0:
        logger.error("No samples given. Exit.")
        sys.exit(-1)

    return (
        sample_dirs,
        sample_names,
        group_filename,
        gff_filename,
        count_filename,
        fastq_filename,
    )


def read_count_info(
    count_filename: Union[str, Path], dirs: List[Union[str, Path]], field_to_use: str
) -> Tuple[str, str]:
    count_info = {}  # key: (sample, PB.1.1) --> count
    count_header = ""
    for name, d in dirs.items():
        with Path(d, count_filename).open() as f:
            while True:
                cur = f.tell()
                line = f.readline().strip()
                if not line.startswith("#"):
                    break
                count_header += line
            f.seek(cur)
            for r in DictReader(f, delimiter="\t"):
                count_info[name, r["pbid"]] = r[field_to_use]
    return count_header, count_info


def chain_split_file(
    ref_gff,
    ref_group,
    ref_name,
    addon_gff,
    addon_group,
    addon_name,
    fuzzy_junction,
    allow_5merge,
    max_3_diff,
    n_chunks,
):
    addon_group_info = sp.MegaPBTree.read_group(addon_group, None)
    recs = []
    tree = OrderedDict()
    i = 0
    for r in GFF.collapseGFFReader(addon_gff):
        if r.chr not in tree:
            tree[r.chr] = {"+": ClusterTree(0, 0), "-": ClusterTree(0, 0)}
        tree[r.chr][r.strand].insert(r.start, r.end, i)
        recs.append(r)
        i += 1

    n = len(recs)
    chunk_size = (n // n_chunks) + (n % n_chunks > 0)

    split_files = []
    i = 0
    counter = 0
    f_gff = open(f"{addon_gff}.split{str(i)}", "w")
    f_group = open(f"{addon_group}.split{str(i)}", "w")
    for v1 in tree.values():
        for strand in ("+", "-"):
            v2 = v1[strand]
            for *_, _indices in v2.getregions():
                for cur in _indices:
                    GFF.write_collapseGFF_format(f_gff, recs[cur])
                    f_group.write(
                        f"{recs[cur].seqid}\t{','.join(addon_group_info[recs[cur].seqid])}\n"
                    )
                    counter += 1
            if counter >= (i + 1) * chunk_size:
                i += 1
                n = f_gff.tell()
                f_gff.close()
                f_group.close()
                if n == 0:  # didn't write any records, delete these
                    Path(f_gff.name).unlink()
                    Path(f_group.name).unlink()
                else:
                    split_files.append((f_gff.name, f_group.name))
                if i >= n_chunks or counter >= len(recs):
                    break
                f_gff = open(f"{addon_gff}.split{str(i)}", "w")
                f_group = open(f"{addon_group}.split{str(i)}", "w")
    if not f_gff.closed:
        n = f_gff.tell()
        f_gff.close()
        f_group.close()
        if n == 0:  # didn't write any records, delete these
            Path(f_gff.name).unlink()
            Path(f_group.name).unlink()
        else:
            split_files.append((f_gff.name, f_group.name))

    result_prefixes = []
    pools = []
    for i, (split_gff, split_group) in enumerate(split_files):
        p = Process(
            target=chain_helper,
            args=(
                ref_gff,
                ref_group,
                split_gff,
                split_group,
                ref_name,
                f"{addon_name}.{str(i)}",
                fuzzy_junction,
                allow_5merge,
                max_3_diff,
            ),
        )
        p.start()
        pools.append(p)
        result_prefixes.append((ref_name, f"{addon_name}.{str(i)}"))
    for p in pools:
        p.join()
    return result_prefixes, split_files


def chain_helper(
    ref_gff: Union[str, Path],
    ref_group: Union[str, Path],
    addon_gff: Union[str, Path],
    addon_group: Union[str, Path],
    name1: str,
    name2: str,
    fuzzy_junction: int,
    allow_5merge: bool,
    max_3_diff: int,
) -> None:
    o = sp.MegaPBTree(
        gff_filename=ref_gff,
        group_filename=ref_group,
        self_prefix=name1,
        internal_fuzzy_max_dist=fuzzy_junction,
        allow_5merge=allow_5merge,
        max_3_diff=max_3_diff,
        fastq_filename=None,
    )
    o.add_sample(
        gff_filename=addon_gff,
        group_filename=addon_group,
        sample_prefix=name2,
        output_prefix=f"tmp_{name2}",
        fastq_filename=None,
    )


def combine_split_chained_results(
    output_prefixes,
    final_prefix,
    ref_gff,
    ref_group,
    ref_name,
    ref_fq,
    addon_gff,
    addon_group,
    addon_name,
    addon_fq,
):
    """
    Each <output_prefix> will have .gff, .group.txt, .mega_info.txt.
    There should be NO overlap between the split files, so clean merge should be possible!

    1. read the .gff files, record the group and mega (id-map) info
    2. sort the total records so can properly put on a unified superPBID
    3. write out the unified result
    4. delete the split files
    """

    # sanity check files are all there
    split_files = []  # tuple of (gff, group, mega)
    for ref_name, o in output_prefixes:
        gff_file = Path(f"tmp_{o}.gff")
        mega_file = Path(f"tmp_{o}.mega_info.txt")
        group_file = Path(f"tmp_{o}.group.txt")
        if not gff_file.exists() or not mega_file.exists() or not group_file.exists():
            raise RuntimeError(
                f"Expects to see {gff_file},{mega_file},{group_file} but one or more files are missing! Abort!"
            )
        split_files.append((ref_name, o, gff_file, group_file, mega_file))

    use_fq = False
    if ref_fq is not None and addon_fq is not None:
        use_fq = True
        ref_fq_dict = {
            r.id.split("|")[0]: r for r in SeqIO.parse(open(ref_fq), "fastq")
        }
        addon_fq_dict = {
            r.id.split("|")[0]: r for r in SeqIO.parse(open(addon_fq), "fastq")
        }

    mega_info = {}  # ref id -> list of matching query_id, or empty list
    split_unmatched = set()

    for (ref_name, split_name, gff_file, group_file, mega_file) in split_files:
        for r in DictReader(open(mega_file), delimiter="\t"):
            if r[ref_name] != "NA":
                if r[ref_name] not in mega_info:
                    mega_info[r[ref_name]] = []
                if r[split_name] != "NA":
                    mega_info[r[ref_name]].append(r[split_name])
            else:  # ref is NA, non-ref is not NA
                split_unmatched.add(r[split_name])

    # make a rec list of matches of (ref_id, addon_id, representative record, combined group info) where rec_ref or ref_addon could be None, but not both
    rec_list = []
    d_ref = {r.seqid: r for r in GFF.collapseGFFReader(ref_gff)}
    d_addon = {r.seqid: r for r in GFF.collapseGFFReader(addon_gff)}

    ref_group_info = sp.MegaPBTree.read_group(ref_group, None)
    addon_group_info = sp.MegaPBTree.read_group(addon_group, None)

    for ref_id, matches in mega_info.items():
        if len(matches) == 0:
            rec_list.append(
                sp.MatchRecord(
                    ref_id=ref_id,
                    addon_id="NA",
                    rec=d_ref[ref_id],
                    members=ref_group_info[ref_id],
                    seqrec=ref_fq_dict[ref_id] if use_fq else None,
                )
            )
        else:
            for addon_id in matches:
                r1 = d_ref[ref_id]
                r2 = d_addon[addon_id]
                if (r1.end - r1.start) > (r2.end - r2.start):
                    rec_list.append(
                        sp.MatchRecord(
                            ref_id=ref_id,
                            addon_id=addon_id,
                            rec=r1,
                            members=ref_group_info[ref_id] + addon_group_info[addon_id],
                            seqrec=ref_fq_dict[ref_id] if use_fq else None,
                        )
                    )
                else:
                    rec_list.append(
                        sp.MatchRecord(
                            ref_id=ref_id,
                            addon_id=addon_id,
                            rec=r2,
                            members=ref_group_info[ref_id] + addon_group_info[addon_id],
                            seqrec=addon_fq_dict[addon_id] if use_fq else None,
                        )
                    )
    for addon_id in split_unmatched:
        rec_list.append(
            sp.MatchRecord(
                ref_id="NA",
                addon_id=addon_id,
                rec=d_addon[addon_id],
                members=addon_group_info[addon_id],
                seqrec=addon_fq_dict[addon_id] if use_fq else None,
            )
        )

    sp.write_reclist_to_gff_n_info(rec_list, final_prefix, ref_name, addon_name, use_fq)
    for (ref_name, split_name, gff_file, group_file, mega_file) in split_files:
        gff_file.unlink()
        group_file.unlink()
        mega_file.unlink()


def chain_samples(
    dirs,
    names,
    group_filename,
    gff_filename,
    count_filename,
    field_to_use="count_fl",
    fuzzy_junction=0,
    allow_5merge=False,
    max_3_diff=100,
    fastq_filename=None,
):
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
        o = sp.MegaPBTree(
            f"tmp_{name}.gff",
            f"tmp_{name}.group.txt",
            self_prefix=f"tmp_{name}",
            internal_fuzzy_max_dist=fuzzy_junction,
            allow_5merge=allow_5merge,
            max_3_diff=max_3_diff,
            fastq_filename=f"tmp_{name}.rep.fq" if fastq_filename is not None else None,
        )
        # chain.append(name) # no need, already done above
    else:  # everything is new, start fresh
        name = names[0]
        d = Path(dirs[name])
        chain = [name]
        o = sp.MegaPBTree(
            d.joinpath(gff_filename),
            d.joinpath(group_filename),
            self_prefix=name,
            internal_fuzzy_max_dist=fuzzy_junction,
            allow_5merge=allow_5merge,
            max_3_diff=max_3_diff,
            fastq_filename=d.joinpath(fastq_filename)
            if fastq_filename is not None
            else None,
        )
        start_i = 1

    for name in names[start_i:]:
        if name.startswith("tmp_"):
            raise AssertionError("trying to add a temp file!")
        d = Path(dirs[name])
        o.add_sample(
            d.joinpath(gff_filename),
            d.joinpath(group_filename),
            sample_prefix=name,
            output_prefix=f"tmp_{name}",
            fastq_filename=d.joinpath(fastq_filename)
            if fastq_filename is not None
            else None,
        )
        o = sp.MegaPBTree(
            f"tmp_{name}.gff",
            f"tmp_{name}.group.txt",
            self_prefix=f"tmp_{name}",
            internal_fuzzy_max_dist=fuzzy_junction,
            allow_5merge=allow_5merge,
            max_3_diff=max_3_diff,
            fastq_filename=f"tmp_{name}.rep.fq" if fastq_filename is not None else None,
        )
        chain.append(name)

    # now recursively chain back by looking at mega_info.txt!!!
    d = {}  # ex: (tmp_1009, PB.1.1) --> mega info dict
    for c in chain[1:]:
        for r in DictReader(open(f"tmp_{c}.mega_info.txt"), delimiter="\t"):
            d[f"tmp_{c}", r["superPBID"]] = r

    with open("all_samples.chained_ids.txt", "w") as f1, open(
        "all_samples.chained_count.txt", "w"
    ) as f2:
        writer1 = DictWriter(f1, fieldnames=["superPBID"] + chain, delimiter="\t")
        writer1.writeheader()

        writer2 = DictWriter(f2, fieldnames=["superPBID"] + chain, delimiter="\t")
        writer2.writeheader()

        reader = DictReader(open(f"tmp_{chain[-1]}.mega_info.txt"), delimiter="\t")
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

            rec1 = {"superPBID": r0["superPBID"]}
            rec2 = {"superPBID": r0["superPBID"]}
            for c in chain:
                rec1[c] = answer[c]
                rec2[c] = str(answer2[c])
            writer1.writerow(rec1)
            writer2.writerow(rec2)

    shutil.copyfile(f"tmp_{chain[-1]}.gff", "all_samples.chained.gff")
    if fastq_filename is not None:
        shutil.copyfile(f"tmp_{chain[-1]}.rep.fq", "all_samples.chained.rep.fq")

    logger.info("Chained output written to:")
    logger.info("all_samples.chained.gff")
    logger.info(f1.name)
    logger.info(f2.name)
    if fastq_filename is not None:
        logger.info("all_samples.chained.rep.fq")


def chain_samples_multithread(
    dirs,
    names,
    group_filename,
    gff_filename,
    count_filename,
    field_to_use="count_fl",
    fuzzy_junction=0,
    allow_5merge=False,
    max_3_diff=100,
    fastq_filename=None,
    cpus=4,
):
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
        first_add = False
    else:  # everything is new, start fresh
        name = names[0]
        chain = [name]
        start_i = 1
        first_add = True

    for addon_name in names[start_i:]:
        assert not addon_name.startswith("tmp_")
        ref_name = chain[-1]
        ref_d = Path(dirs[ref_name])
        if first_add:
            ref_gff = ref_d.joinpath(gff_filename)
            ref_group = ref_d.joinpath(group_filename)
            ref_fq = (
                ref_d.joinpath(fastq_filename) if fastq_filename is not None else None
            )
        else:
            ref_name = f"tmp_{ref_name}"
            ref_gff = f"{ref_name}.gff"
            ref_group = f"{ref_name}.group.txt"
            ref_fq = f"{ref_name}.rep.fq" if fastq_filename is not None else None
        addon_d = Path(dirs[addon_name])
        addon_gff = addon_d.joinpath(gff_filename)
        addon_group = addon_d.joinpath(group_filename)
        addon_fq = (
            addon_d.joinpath(fastq_filename) if fastq_filename is not None else None
        )
        split_outs, split_ins = chain_split_file(
            ref_gff=ref_gff,
            ref_group=ref_group,
            ref_name=ref_name,
            addon_gff=addon_gff,
            addon_group=addon_group,
            addon_name=addon_name,
            fuzzy_junction=fuzzy_junction,
            allow_5merge=allow_5merge,
            max_3_diff=max_3_diff,
            n_chunks=cpus,
        )

        combine_split_chained_results(
            split_outs,
            final_prefix=f"tmp_{addon_name}",
            ref_gff=ref_gff,
            ref_group=ref_group,
            ref_name=ref_name,
            ref_fq=ref_fq,
            addon_gff=addon_gff,
            addon_group=addon_group,
            addon_name=addon_name,
            addon_fq=addon_fq,
        )

        chain.append(addon_name)
        for in_gff_split, in_group_split in split_ins:
            Path(in_gff_split).unlink()  # remove the split gff
            Path(in_group_split).unlink()

        first_add = False

    # now recursively chain back by looking at mega_info.txt!!!
    d = {}  # ex: (tmp_sample1, PB.1.1) --> mega info dict
    for c in chain[1:]:
        for r in DictReader(open(f"tmp_{c}.mega_info.txt"), delimiter="\t"):
            d[f"tmp_{c}", r["superPBID"]] = r

    with open("all_samples.chained_ids.txt", "w") as f1, open(
        "all_samples.chained_count.txt", "w"
    ) as f2:
        writer1 = DictWriter(f1, fieldnames=["superPBID"] + chain, delimiter="\t")
        writer1.writeheader()

        writer2 = DictWriter(f2, fieldnames=["superPBID"] + chain, delimiter="\t")
        writer2.writeheader()

        reader = DictReader(open(f"tmp_{chain[-1]}.mega_info.txt"), delimiter="\t")
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

            rec1 = {"superPBID": r0["superPBID"]}
            rec2 = {"superPBID": r0["superPBID"]}
            for c in chain:
                rec1[c] = answer[c]
                rec2[c] = str(answer2[c])
            writer1.writerow(rec1)
            writer2.writerow(rec2)

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
        help="Which count field to use for chained sample (default: count_fl)",
    ),
    fuzzy_junction: int = typer.Option(
        0,
        show_default=False,
        help="Max allowed distance in junction to be considered identical (default: 0 bp)",
    ),
    allow_5merge: bool = typer.Option(
        True,
        "--dun-merge-5-shorter",
        show_default=False,
        help="Don't collapse shorter 5' transcripts (default: off)",
    ),  # store_false
    max_3_diff: int = typer.Option(
        30, show_default=False, help="Maximum 3' difference allowed (default: 30bp)"
    ),
    cpus: int = typer.Option(
        8,
        show_default=False,
        help="Number of CPUs to use for multi-threading (default: 8)",
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
    chain_samples_multithread(
        sample_dirs,
        sample_names,
        group_filename,
        gff_filename,
        count_filename,
        field_to_use.value,
        fuzzy_junction,
        allow_5merge,
        max_3_diff,
        fastq_filename,
        cpus=cpus,
    )


if __name__ == "__main__":
    typer.run(main)
