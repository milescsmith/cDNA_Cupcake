#!/usr/bin/env python
__author__ = "etseng@pacb.com"

import sys
from csv import DictWriter
from enum import Enum
from pathlib import Path

import typer
from Bio import SeqIO
from bx.intervals import Interval, IntervalTree
from cupcake.sequence import BED, BioReaders
from cupcake.sequence.GFF import collapseGFFReader


class index_base(int, Enum):
    zero = (0,)
    one = 1


app = typer.Typer(
    name="cupcake.targeted.calc_probe_hit_from_sam",
    help="Calculate Probe Hit from SAM alignment + probe BED",
)


def get_probe_hit(tree, gene_info, r, is_gtf=False):
    """
    Given a dict tree (from read_probe_bed) and a GMAP SAM record
    Go through each exon and find probes that hit it

    Return: (number of probes hit), (total number of bases overlapping with probes), (genes seen)
    """
    probes_seen = set()
    genes_seen = set()
    base_hit = 0
    if is_gtf:
        r.sID, r.segments = r.chr, r.ref_exons
    if r.sID not in tree:
        return 0, 0, set()
    for e in r.segments:
        hits = tree[r.sID].find(e.start, e.end)
        if len(hits) == 0:
            continue
        for i, strand, intl in hits:
            if (strand is None) or (strand == r.strand):
                probes_seen.add(i)
                genes_seen.add(gene_info[i])
                base_hit += min(e.end, intl.end) - max(e.start, intl.start)
    return len(probes_seen), base_hit, genes_seen


def read_probe_bed(bed_filename, start_base=0, end_base=1):
    """
    Read a probe BED file <chrom>, <start>, <end>, [optional score], [optional strand]
    Return dict of chrom --> IntervalTree w/ data=(index, interval)
    """
    tree = {}
    gene_info = {}
    i = 0
    reader = BED.SimpleBEDReader(bed_filename, start_base, end_base)
    for r in reader:
        if r.chr not in tree:
            tree[r.chr] = IntervalTree()
        tree[r.chr].add(r.start, r.end, (i, r.strand, Interval(r.start, r.end)))
        if r.name is not None:
            gene_info[i] = r.name
        i += 1
    return tree, gene_info


def calc_ontarget_rate(
    tree, gene_info, input_fasta, is_gtf, sam_or_gtf, output_filename=None
):

    feature = (
        "fasta"
        if input_fasta.upper().endswith(".FA") or input_fasta.upper().endswith(".FASTA")
        else "fastq"
    )
    query_len_dict = {r.id: len(r.seq) for r in SeqIO.parse(open(input_fasta), feature)}

    if output_filename is None:
        f = sys.stdout
    else:
        f = Path(output_filename)

    FIELDS = ["read_id", "read_len", "num_probe", "num_base_overlap", "loci", "genes"]
    writer = DictWriter(f, FIELDS, delimiter="\t")
    writer.writeheader()

    if is_gtf:
        reader = collapseGFFReader(sam_or_gtf)
        for r in reader:
            num_probe, base_hit, genes_seen = get_probe_hit(tree, gene_info, r, is_gtf)
            rec = {
                "read_id": r.seqid,
                "read_len": "NA",
                "num_probe": num_probe,
                "num_base_overlap": base_hit,
                "loci": f"{r.chr}:{r.start}-{r.end}",
                "genes": ",".join(genes_seen),
            }
            writer.writerow(rec)
    else:
        reader = BioReaders.GMAPSAMReader(
            sam_or_gtf, True, query_len_dict=query_len_dict
        )
        for r in reader:
            if r.sID == "*":
                continue
            num_probe, base_hit, genes_seen = get_probe_hit(tree, gene_info, r, is_gtf)
            rec = {
                "read_id": r.qID,
                "read_len": r.qLen,
                "num_probe": num_probe,
                "num_base_overlap": base_hit,
                "loci": f"{r.sID}:{r.sStart}-{r.sEnd}",
                "genes": ",".join(genes_seen),
            }
            writer.writerow(rec)


@app.command()
def main(
    bed_filename: str = typer.Argument(...),
    input_fasta_or_fastq: str = typer.Argument(...),
    sam_or_gtf: str = typer.Argument(...),
    gtf: bool = typer.Option(False, help="Input is GTF instead of SAM"),
    start_base: index_base = typer.Option(..., help="Start base is 0 or 1-based index"),
    end_base: index_base = typer.Option(..., help="End base is 0 or 1-based index"),
    output: str = typer.Option(
        "stdout", "-o", "--output", help="Output filename (default: stdout)"
    ),
):

    tree, gene_info = read_probe_bed(bed_filename, int(start_base), int(end_base))

    calc_ontarget_rate(
        tree,
        gene_info,
        input_fasta_or_fastq,
        gtf,
        sam_or_gtf,
        output,
    )


if __name__ == "__main__":
    typer.run(main)
