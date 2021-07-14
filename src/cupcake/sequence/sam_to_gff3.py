#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Convert a SAM file into GFF3 format (https://uswest.ensembl.org/info/website/upload/gff3.html)
 that also includes additional information that GMAP GFF3 for the 'mRNA' type:
   --- coverage
   --- identity
   --- matches
   --- mismatches
   --- indels


ex: from GMAP

6   cow_hereford    mRNA    109018744   109018861   .   +   .   \
     ID=myid.mrna1;Name=myid;Parent=myid.path1;coverage=29.5;identity=98.3;matches=116;mismatches=2;indels=0;unknowns=0
"""


from collections import Counter
from pathlib import Path
from typing import Optional

import typer

# from gtfparse.write_gtf import df_to_gtf
from BCBio import GFF as BCBio_GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from cupcake import version_callback
from cupcake.logger import cupcake_logger as logger
from cupcake.sequence.BioReaders import GMAPSAMReader

app = typer.Typer(
    name="cupcake.sequence.sam_to_gff3",
    help="Convert SAM to GFF3 format using BCBio GFF",
)


def convert_sam_rec_to_gff3_rec(r, source, qid_index_dict=None):
    """
    :param r: GMAPSAMRecord record
    :param qid_seen: list of qIDs processed so far -- if redundant, we have to put a unique suffix
    :return SeqRecord ready to be written as GFF3
    """
    if r.sID == "*":
        logger.info(f"Skipping {r.qID} because unmapped.")
        return None
    t_len = sum(e.end - e.start for e in r.segments)
    seq = Seq("A" * t_len)  # DO NOT CARE since sequence is not written in GFF3
    rec = SeqRecord(seq, r.sID)
    strand = 1 if r.flag.strand == "+" else -1

    # indels = r.num_ins + r.num_del
    # mismatches = r.num_nonmatches
    # matches = r.num_mat_or_sub - r.num_nonmatches

    if qid_index_dict is not None:
        if r.qID in qid_index_dict:
            qid_index_dict[r.qID] += 1
            r.qID += f"_dup{str(qid_index_dict[r.qID])}"
        else:
            qid_index_dict[r.qID] += 1

    gene_qualifiers = {"source": source, "ID": r.qID, "Name": r.qID}  # for gene record
    #    mRNA_qualifiers = {"source": source, "ID": r.qID+'.mRNA', "Name": r.qID+'.mRNA', "Parent": r.qID,
    #                       "coverage": "{0:.2f}".format(r.qCoverage*10**2) if r.qCoverage is not None else "NA",
    #                       "identity": "{0:.2f}".format(r.identity*10**2),
    #                       "matches": matches, "mismatches": mismatches, "indels": indels}

    # gene line, one per record
    top_feature = SeqFeature(
        FeatureLocation(r.sStart, r.sEnd),
        type="gene",
        strand=strand,
        qualifiers=gene_qualifiers,
    )
    # mRNA line, one per record
    top_feature.sub_features = (
        []
    )  # top_feature.sub_features = [SeqFeature(FeatureLocation(r.sStart, r.sEnd), type="mRNA", strand=strand, qualifiers=mRNA_qualifiers)]

    # exon lines, as many exons per record
    for i, e in enumerate(r.segments):
        _id = f"{r.qID}.exon{i+1}"
        exon_qual = {"source": source, "ID": _id, "Name": _id}
        top_feature.sub_features.append(
            SeqFeature(
                FeatureLocation(e.start, e.end),
                type="exon",
                strand=strand,
                qualifiers=exon_qual,
            )
        )
    rec.features = [top_feature]
    return rec


def convert_sam_to_gff3(sam_filename: str, output_gff3: str, source, q_dict=None):
    qid_index_dict = Counter()
    with open(output_gff3, "w") as f:
        recs = [
            convert_sam_rec_to_gff3_rec(r0, source, qid_index_dict)
            for r0 in GMAPSAMReader(sam_filename, True, query_len_dict=q_dict)
        ]
        BCBio_GFF.write([x for x in recs if x is not None], f)


@app.command(name="")
def main(
    sam_filename: str = typer.Argument(...),
    input_fasta: Optional[str] = typer.Option(
        None,
        "--input_fasta",
        "-i",
        help="(Optional) input fasta. If given, coverage will be calculated.",
    ),
    source: str = typer.Option(
        ..., "--source", "-s", help="source name (ex: hg38, mm10)"
    ),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
):
    sam_filename = Path(sam_filename)

    if sam_filename.suffix != (".sam"):
        raise RuntimeError("Only accepts files ending in .sam. Abort!")

    prefix = sam_filename.stem
    output_gff3 = f"{prefix}.gff3"

    q_dict = None
    if input_fasta is not None:
        q_dict = {r.id: len(r.seq) for r in SeqIO.parse(open(input_fasta), "fasta")}

    with open(output_gff3, "w") as f:
        recs = [
            convert_sam_rec_to_gff3_rec(r0, source)
            for r0 in GMAPSAMReader(sam_filename, True, query_len_dict=q_dict)
        ]
        BCBio_GFF.write([x for x in recs if x is not None], f)

    logger.info(f"Output written to {output_gff3}.")


if __name__ == "__main__":
    typer.run(main)
