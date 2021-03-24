#!/usr/bin/env python
"""
A temporary CSV file for isoseq3 (v3.4+) dedup output

INPUT: dedup.fasta
OUTPUT: dedup.info.csv
"""
import re

from Bio import SeqIO
from Bio.Seq import Seq

rex = re.compile(r"(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+);XC=(\S+)")
rex_umi_only = re.compile(r"(\S+) full_length_coverage=(\d+);length=(\d+);XM=(\S+)")

reader = SeqIO.parse(open("dedup.fasta"), "fasta")
with open("dedup.info.csv", "w") as f:
    f.write("id\tUMI\tUMIrev\tBC\tBCrev\tlength\tcount\n")
    for r in reader:
        m = rex.match(r.description)
        if m is not None:
            _id, _count, _len, _umi, _bc = m.groups()
            f.write(
                f"{_id}\t{_umi}\t{Seq(_umi).reverse_complement()}\t{_bc}\t{Seq(_bc).reverse_complement()}\t{_len}\t{_count}\n"
            )
        else:
            m = rex_umi_only.match(r.description)
            _id, _count, _len, _umi = m.groups()
            _bc = "NA"
            f.write(
                f"{_id}\t{_umi}\t{Seq(_umi).reverse_complement()}\t{_bc}\t{Seq(_bc).reverse_complement()}\t{_len}\t{_count}\n"
            )
