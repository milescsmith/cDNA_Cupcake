#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Should be used *after* summarize_sample_GFF_junctions.py is called to generate junction report.

Looks through the junction reports and scrub it, retaining only junctions that
 meet one or more of the following criteria:

(a). be the only junction in that "label" <-- i.e. you must keep this junction because there's no other similar ones to it
(b). have previous annotation (annotation='Y') <-- always trust reference annotation
(c). have at least X sample supporting it (`num_sample>=T` or `num_transcript>=S`) where `S`, `T` is user-defined
(d). additionally, --accept-all-canonical <-- meaning all canonical junctions are accepted

Output: scrubbed.junctions.bed
"""
import os
import sys
from collections import defaultdict
from csv import DictReader, DictWriter
from pathlib import Path
from typing import Union, Optional, Tuple, List, Dict, Any
import typer
from Bio import SeqIO
from bx.intervals import Interval, IntervalTree

import cupcake.sequence.GFF as GFF
import cupcake.tofu.counting.chain_samples as sp
from cupcake.logging import cupcake_logger as logger

fields_to_add = [
    "count_fl",
    "count_nfl",
    "count_nfl_amb",
    "norm_fl",
    "norm_nfl",
    "norm_nfl_amb",
]


app = typer.Typer(name="cupcake.tofu.counting.scrub_sample_GFF_junctions")


def read_junction_report(filename: Union[str, Path]) -> Tuple[str, List[str]]:
    """
    tab-delimited with header:
        chr     left    right   strand  num_transcript  num_sample      genome  annotation      label

    return: dict of label --> records with that label
    """
    with open(filename, "r") as fi:
        reader = DictReader(open(fi), delimiter="\t")
        r = next(reader)

        cur_label, cur = r["label"], [r]
        for r in reader:
            if r["label"] != cur_label:
                yield cur_label, cur
                cur_label, cur = r["label"], [r]
            else:
                cur.append(r)
        yield cur_label, cur


def scrub_junction_by_label(
    junctions: List[Dict[str, Any]], min_sample: int = 2, min_transcript: int = 2, accept_all_canonical: bool = False
) -> List[Dict[str, Any]]:
    """
    input: a list of junctions (records) that are of the same label
    output: a list of "scrubbed" junctions
    """

    def pass_criteria(r: Dict[str, Any]) -> bool:
        if int(r["num_sample"]) >= min_sample:
            return True
        if int(r["num_transcript"]) >= min_transcript:
            return True
        if accept_all_canonical and r["genome"] == "GT-AG":
            return True
        if r["annotation"] == "Y":
            return True

    if len(junctions) == 1:
        # only one junction in this label, must accept!
        return junctions

    good = []
    for r in junctions:
        if pass_criteria(r):
            good.append(r)

    if (
        len(good) == 0
    ):  # nothing passed criteria! must pick one. we will pick one that has the highest transcript support
        junctions.sort(key=lambda r: int(r["num_transcript"]), reverse=True)
        good.append(junctions[0])
    return good


def find_best_match_junction(
    tree: IntervalTree,
    donor: int,
    accep: int,
    max_diff: int = 20,
) -> Optional[Interval]:
    """
    donor, accept -- both should be 0-based
    """
    hits = tree.find(donor, accep)
    if len(hits) == 0:
        return None
    elif len(hits) == 1:
        if hits[0].start - donor > max_diff or hits[0].end - accep > max_diff:
            return None
        return hits[0]
    else:  # multiple hits, find the closest one
        diff = []
        for h in hits:
            if h.start - donor > max_diff or h.end - accep > max_diff:
                continue
            diff.append((abs(h.start - donor) + abs(h.end - accep), h))
        diff.sort(key=lambda x: x[0])
        return diff[0][1]


def scrub_ref_exons(r: Dict[str, Any], tree: IntervalTree) -> Optional[List[Interval]]:
    n = len(r.ref_exons)
    new_ref_exons = []
    cur_start = r.ref_exons[0].start
    for i in range(n - 1):
        donor = r.ref_exons[i].end - 1  # make it 0-based
        accep = r.ref_exons[i + 1].start  # start is already 0-based
        match = find_best_match_junction(tree[r.chr, r.strand], donor, accep)
        if match is None:
            logger.info(f"donor-acceptor site {r.chr},{r.strand},{donor}-{accep} has no hit in tree!")
            return None

        new_ref_exons.append(Interval(cur_start, match.start + 1))
        cur_start = match.end
    new_ref_exons.append(Interval(cur_start, r.ref_exons[-1].end))
    return new_ref_exons


def read_scrubbed_junction_to_tree(junction_filename: Union[str, Path]) -> IntervalTree:
    tree = defaultdict(lambda: IntervalTree())
    with open(junction_filename) as f:
        if not f.readline().startswith("track"):
            f.seek(0)
        for line in f:
            raw = line.strip().split("\t")
            if len(raw) == 4:
                chrom, left, right, strand = raw
            elif len(raw) == 6:
                chrom, left, right, _name, _count, strand = raw
            else:
                raise Exception(f"Expects junction BED file to have either 4 or 6 columns! Saw {len(raw)}!")
            left, right = int(left), int(right)  # already 0-based start, 0-based end
            tree[chrom, strand].add(left, right, Interval(left, right))
    return tree


def scrub_junctions(
    report_filename: Union[str, Path],
    output_filename: Union[str, Path],
    min_sample: int,
    min_transcript: int,
    accept_all_canonical: bool
) -> IntervalTree:
    tree = defaultdict(lambda: IntervalTree())
    with open(output_filename, "w") as f:
        for _, junctions in read_junction_report(report_filename):
            good = scrub_junction_by_label(
                junctions, min_sample, min_transcript, accept_all_canonical
            )
            for r in good:
                a, b = int(r["left"]), int(r["right"])  # 0-based start, 0-basde end
                f.write(
                    f"{r['chr']}\t{r['left']}\t{r['right']}\t{r['strand']}\n"
                )
                tree[r["chr"], r["strand"]].add(a, b, Interval(a, b))
    return tree


def scrub_sample_GFFs(
    sample_dirs: Dict[str, str],
    gff_filename: Union[str, Path],
    count_filename: Union[str, Path],
    group_filename: Union[str, Path],
    fastq_filename: Union[str, Path],
    output_prefix: str,
    tree: IntervalTree,
) -> None:

    for _, d in sample_dirs.items():
        with Path(d, f"{output_prefix}.gff.tmp").open("w") as outf:
            for r in GFF.collapseGFFReader(os.path.join(d, gff_filename)):
                n = len(r.ref_exons)
                if n == 1:
                    GFF.write_collapseGFF_format(outf, r)

                new_ref_exons = scrub_ref_exons(r, tree)
                if new_ref_exons is None:
                    print("No changes made due to error:", r.seqid, file=sys.stderr)
                else:
                    # print "before:", r.ref_exons
                    # print "after :", new_ref_exons
                    r.ref_exons = new_ref_exons
                GFF.write_collapseGFF_format(outf, r)
        cleanup_scrubbed_files_redundancy(
            outf.name,
            Path(d, group_filename),
            Path(d, count_filename),
            Path(d, fastq_filename) if fastq_filename is not None else None,
            Path(d, output_prefix),
        )


def read_count_file(count_filename: Union[str, Path]) -> Tuple[Dict[str, Dict[str, Any]]]:
    with open(count_filename) as f:
        count_header = ""
        while True:
            cur_pos = f.tell()
            line = f.readline()
            if not line.startswith("#"):
                f.seek(cur_pos)
                break
            else:
                count_header += line
        d = {r["pbid"]: r for r in DictReader(f, delimiter="\t")}
    return d, count_header


def read_group_file(group_filename):
    group_info = {}  # key: PB.1.1 --> list of HQ isoforms
    for line in open(group_filename):
        a, b = line.strip().split("\t")
        group_info[a] = b
    return group_info


def cleanup_scrubbed_files_redundancy(
    gff_filename: Union[str, Path],
    group_filename: Union[str, Path],
    count_filename: Union[str, Path],
    fastq_filename: Union[str, Path],
    output_prefix: str
):

    junction_seen = defaultdict(
        lambda: defaultdict(lambda: [])
    )  # key (chr,strand) -> dict of (series of junctions) -> record
    for r in GFF.collapseGFFReader(gff_filename):
        n = len(r.ref_exons)
        if n == 1:
            junc_str = f"{str(r.start)},{str(r.end)}"
            junction_seen[r.chr, r.strand][junc_str] = [r]
        else:
            junc_str = ",".join(
                f"{str(r.ref_exons[i].end)},{str(r.ref_exons[i + 1].start)}"
                for i in range(n - 1)
            )
            junction_seen[r.chr, r.strand][junc_str].append(r)

    # write out cleaned GFF
    with open(f"{output_prefix}.gff", "w") as outf, open(f"{output_prefix}.merged_ids.txt", "w") as outf2:
        merged = {}
        keys = list(junction_seen.keys())
        keys.sort()
        for k in keys:
            for bunch in junction_seen[k].values():
                if len(bunch) == 1:  # just one record, write it out
                    r = bunch[0]
                    GFF.write_collapseGFF_format(outf, r)
                    merged[r.seqid] = [r.seqid]
                else:
                    # find the representative
                    r = bunch[0]
                    for r2 in bunch[1:]:
                        if r2.end - r2.start > r.end - r.start:
                            r = r2
                    GFF.write_collapseGFF_format(outf, r)
                    merged[r.seqid] = [x.seqid for x in bunch]
                outf2.write(f"{r.seqid}\t{','.join(merged[r.seqid])}\n")

    count_d, count_header = read_count_file(count_filename)
    # write out count file
    with open(f"{output_prefix}.abundance.txt", "w") as outf:
        outf.write(count_header)
        writer = DictWriter(
            outf,
            fieldnames=[
                "pbid",
                "count_fl",
                "count_nfl",
                "count_nfl_amb",
                "norm_fl",
                "norm_nfl",
                "norm_nfl_amb",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for pbid, bunch in merged.items():
            # combine the counts
            r = count_d[bunch[0]]
            r["pbid"] = pbid
            for field in fields_to_add:
                r[field] = float(r[field])
            for _id in bunch[1:]:
                for field in fields_to_add:
                    r[field] += float(count_d[_id][field])
            writer.writerow(r)

    group_info = read_group_file(group_filename)
    # write out group file
    with open(f"{output_prefix}.group.txt", "w") as outf:
        for pbid, bunch in merged.items():
            # combine the groups
            g = [group_info[bunch[0]]]
            for _id in bunch[1:]:
                g.append(group_info[_id])
            outf.write(f"{pbid}\t{','.join(g)}\n")

    # write out fastq file if present
    if fastq_filename is not None:
        with open(f"{output_prefix}.rep.fq", "w") as outf:
            for r in SeqIO.parse(open(fastq_filename), "fastq"):
                if r.id.split("|")[0] in merged or r.id in merged:
                    SeqIO.write(r, outf, "fastq")

    logger.info(f"scrubbed files written: {output_prefix}.gff, {output_prefix}.group.txt, {output_prefix}.abundance.txt, {output_prefix}.merged_ids.txt")


@app.command(name="")
def main(
    sample_config: str = typer.Argument(...),
    summary_report: str = typer.Argument(...),
    output_prefix: str = typer.Argument(...),
    min_sample: int = typer.Option(1, "-S", help="Minimum number of samples as evidence (default: 1)"),
    min_transcript: int = typer.Option(2, "-T", help="Minimum number of transcripts as evidence (default: 2)"),
    # parser.add_argument("-C", "--accept_all_canonical", action="store_true", default=False, help="Accept all canonical jucntions (default: false)")
    scrubbed_junction_file: Union[str, Path] = typer.Option(help="Scrubbed junction bed --- if given, directly use it to scrub GFFs.")
):
    (
        sample_dirs,
        sample_names,
        group_filename,
        gff_filename,
        count_filename,
        fastq_filename,
    ) = sp.read_config(sample_config)

    report_filename = summary_report

    if scrubbed_junction_file is None:
        output_filename = f"{output_prefix}.scrubbed.junction.bed"
        tree = scrub_junctions(
            report_filename, output_filename, min_sample, min_transcript, True
        )
        logger.info(f"Scrubbed junction written to: {output_filename}")
    else:
        output_filename = scrubbed_junction_file
        logger.info(f"Reading scrubbed junction file: {output_filename}")
        tree = read_scrubbed_junction_to_tree(output_filename)

    scrub_sample_GFFs(
        sample_dirs,
        gff_filename,
        count_filename,
        group_filename,
        fastq_filename,
        output_prefix,
        tree,
    )


if __name__ == "__main__":
    typer.app(main)
