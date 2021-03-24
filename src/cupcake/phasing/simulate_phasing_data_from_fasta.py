#!/usr/bin/env python
__author__ = "etseng@pacb.com"

"""
Input: fasta, each sequence is an individual gene

given an input transcript sequence ---
-> randomly select N positions to have a SNP, N = <length> / 300
-> generate two haplotypes, one is the original, one is mutated
(store this in Haplotypes)
-> for each haplotype, simulate 100X sequences at 0%, 1%, 2%, 3% (each is a diff output fasta)

"""
import random
from collections import Counter
from pathlib import Path
from typing import List

import typer
from Bio import SeqIO
from cupcake.logging import cupcake_logger as logger
from cupcake.phasing.io.VariantPhaser import Haplotypes
from cupcake.simulate.simulate import sim_seq

SNP_INTERVAL = 300  # 1 snp per 300 bp
base_choices = {
    "A": ("T", "C", "G"),
    "T": ("A", "C", "G"),
    "C": ("A", "T", "G"),
    "G": ("A", "T", "C"),
}


app = typer.Typer(name="simulate_phasing_data_from_fasta")


def simulate_phasing_data(
    seq0: str,
    err_sub: float,
    ploidity: int = 2,
    copies=List[int],
    write_fastq: bool = False,
    working_dir: Path = Path().cwd(),
) -> None:
    """
    :param seq0: transcript sequence
    :param err_sub: error prob
    :param ploidity: ploidity (how many alleles)
    :param copies: list of how many copies to simulate per allele   ex: [10, 10] means 10 for each allele
    """
    n = len(seq0)
    var_positions = random.sample(list(range(n)), n / SNP_INTERVAL)
    var_positions.sort()

    new_seqs = [seq0]
    for ignore in range(ploidity - 1):
        num_tries = 0
        while (
            num_tries < 10
        ):  # after 10 attemps, give up simulating enough divergent seqs (this happens for very short haps)
            num_tries += 1
            new_seq = list(seq0)
            for p in var_positions:
                new_seq[p] = random.choice(base_choices[seq0[p]])
            new_seq = "".join(new_seq)
            if new_seq not in new_seqs:
                new_seqs.append("".join(new_seq))
                break

    count_of_vars_by_pos = {}
    for p in var_positions:
        count_of_vars_by_pos[p] = Counter()  # base --> count of base at pos
    for i, seq in enumerate(new_seqs):
        for p in var_positions:
            count_of_vars_by_pos[p][seq[p]] += copies[i]

    # 1. write out fake.fasta
    working_dir.joinpath("fake.fasta").write_text(f">fake\n{seq0}\n")

    # 2. write fake.mapping.txt
    for i in range(n):
        working_dir.joinpath("fake.mapping.txt").write_text(f"{i},fake,{i}\n")

    # simulate CCS reads
    if write_fastq:
        f = working_dir.joinpath("ccs.fastq")
    else:
        f = working_dir.joinpath("ccs.fasta")
    working_dir.joinpath("fake.read_stat.txt").write_text(
        "id\tlength\tis_fl\tstat\tpbid\n"
    )
    profile = [err_sub, err_sub, err_sub, 1.0]
    for i, copy in enumerate(copies):
        if i >= len(new_seqs):
            break
        for c in range(copy):
            cur_seq, cur_qv = sim_seq(new_seqs[i], profile)
            cur_id = f"hap{i + 1}_{c}"
            if write_fastq:
                f.write_text(f"@{cur_id}\n{cur_seq}\n")
                f.write_text(f"+\n{cur_qv}\n")
            else:
                f.write_text(f">{cur_id}\n{cur_seq}\n")
            working_dir.joinpath("fake.read_stat.txt").write_text(
                f"{cur_id}\t{n}\tY\tunique\tPB.0.0\n"
            )

    ref_at_pos = {p: seq0[p] for p in var_positions}
    hap_obj = Haplotypes(var_positions, ref_at_pos, count_of_vars_by_pos)
    for seq in new_seqs:
        hap = "".join(seq[p] for p in var_positions)
        hap_obj.match_or_add_haplotype(hap)

    hap_obj.get_haplotype_vcf_assignment()
    hap_obj.write_haplotype_to_vcf("fake.mapping.txt", {}, "fake.answer")


@app.command(name="")
def main(
    fasta_filename: str = typer.Arguement(
        ..., help="Fasta file from which to simulate phasing data from."
    ),
    ploidity: int = typer.Option(2, "-p"),
    err_sub: float = typer.Option(...),
    copies: str = typer.Option(...),
    write_fastq: bool = typer.Option(False),
) -> None:

    assert 2 <= ploidity <= 6

    copies = list(map(int, copies.split(",")))
    assert len(copies) == ploidity

    for r in SeqIO.parse(open(fasta_filename), "fasta"):
        d2 = r.id.split("|")[0]
        logger.info(f"making {d2}")
        Path(d2).mkdir(parents=True, exist_ok=True)
        simulate_phasing_data(
            seq0=r.seq.tostring(),
            err_sub=err_sub,
            ploidity=ploidity,
            copies=copies,
            write_fastq=write_fastq,
            working_dir=d2,
        )


if __name__ == "__main__":
    typer.run(main)
