__author__ = "etseng@pacb.com"
from pathlib import Path

from Bio import SeqIO
from cupcake.utils import OpenFile


def check_ids_unique(fa_or_fq_filename: str, is_fq: bool = False) -> None:
    """
    Confirm that a FASTA/FASTQ file has all unique IDs
    (used probably by collapse or fusion finding script)
    """
    seen = set()
    fa_or_fq_filename = Path(fa_or_fq_filename)

    for r in SeqIO.parse(OpenFile(fa_or_fq_filename), "fastq" if is_fq else "fasta"):
        if r.id in seen:
            raise Exception(f"Duplicate id {r.id} detected. Abort!")
        seen.add(r.id)
