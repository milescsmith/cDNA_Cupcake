#!/usr/bin/env python
import random
import sys
from collections import defaultdict
from pathlib import Path

import typer
from Bio import SeqIO

from cupcake.logging import cupcake_logger as logger

simType = ["sub", "ins", "del", "match"]
simTypeSize = 4


app = typer.Typer(
    name="cupcake.simulate.simulate",
    help="Simple error simulation",
)


def throwdice(profile):
    dice = random.random()
    for i in range(simTypeSize):
        if dice < profile[i]:
            return simType[i]


def sim_start(ntimes, profile):
    start = defaultdict(lambda: 0)
    acc = defaultdict(lambda: 0)
    for _ in range(ntimes):
        curpos = 0
        while True:
            simtype = throwdice(profile)
            acc[type] += 1
            if simtype == "match":
                start[curpos] += 1
                break
            elif simtype == "ins":
                # insertion occurred, with 1/4 chance it'll match
                if random.random() < 0.25:
                    start[curpos] += 1
                    break
            # if is sub or del, just advance cursor forward
            curpos += 1
    return start, acc


def sim_seq(seq, profile):
    """
    :param seq: sequence to simulate from
    :param profile: accumulative prob vector for simType, ex: [0.01, 0.05, 0.07, 1.] means 1% sub, 4% ins, 2% del
    :return: simulated sequence, qv string ('!' for err, ']' for no err, currently don't handle deletion)
    """
    nucl = {"A", "T", "C", "G"}
    sim = ""
    qv = ""  # qv string,
    for _, s in enumerate(seq):
        while True:
            simtype = throwdice(profile)
            if simtype == "match":
                sim += s
                qv += "]"
                break
            elif simtype == "ins":
                # insertion occurred, with 1/4 chance it'll match
                choice = random.sample(nucl, 1)[0]
                sim += choice
                qv += "!"
            elif simtype == "sub":  # anything but the right one
                choice = random.sample(nucl.difference([s]), 1)[0]
                sim += choice
                qv += "!"
                break
            elif simtype == "del":  # skip over this
                break
            else:
                raise KeyError(f"Invalid type {simtype}")

    return sim, qv


@app.command(name="")
def main(
    fasta_filename: str = typer.Argument(...),
    output_prefix: str = typer.Argument(...),
    copy: int = typer.Option(
        1,
        help="Number of copies to simulate per input sequence (default: 1)",
    ),
    ins: float = typer.Option(
        0,
        "--ins",
        "-i",
        help="Insert error rate [0-1] (default: 0)",
    ),
    dele: float = typer.Option(
        0,
        "--dele",
        "-d",
        help="Deletion error rate [0-1] (default: 0)",
    ),
    sub: float = typer.Option(
        0,
        "--sub",
        "-s",
        help="Substitution error rate [0-1] (default: 0)",
    ),
) -> None:
    if sub < 0 or sub > 1:
        logger.error("Substitution error must be between 0-1!")
        sys.exit(-1)
    if ins < 0 or ins > 1:
        logger.error("Insertion error must be between 0-1!")
        sys.exit(-1)
    if dele < 0 or dele > 1:
        logger.error("Deletion error must be between 0-1!")
        sys.exit(-1)

    if sub + ins + dele > 1:
        logger.error("Total sub+ins+del error cannot exceed 1!")
        sys.exit(-1)

    profile = [sub, sub + ins, ins + dele, 1.0]

    fasta_filename = Path(fasta_filename)
    idpre = output_prefix

    ith = 0
    for r in SeqIO.parse(open(fasta_filename), "fasta"):
        for _ in range(copy):
            ith += 1
            print(
                f">{idpre}_{ith}_{r.id[:r.id.find('|')]}\n{sim_seq(r.seq.tostring(), profile)}"
            )


if __name__ == "__main__":
    typer.run(main)
