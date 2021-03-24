#!/usr/bin/env python
import math
import random
from collections import defaultdict
from csv import DictReader
from pathlib import Path
from typing import List, Optional, Tuple, Union

import typer
from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.annotation.subsample")


def get_counts(
    count_filename: Union[str, Path],
    min_fl_count=2,
    key="id",
    min_len=None,
    max_len=None,
) -> Tuple[int, List[int]]:
    total = 0
    count_d = defaultdict(lambda: 0)
    for r in DictReader(open(count_filename), delimiter="\t"):
        _len = int(r["length"])
        if min_len is not None and _len < min_len:
            continue
        if max_len is not None and _len > max_len:
            continue
        c = int(r["fl_count"])
        if c >= min_fl_count:
            count_d[r[key]] += c
            total += c

    counts = []
    for k, v in count_d.items():
        counts += [k] * v

    return total, counts


def subsample(
    total: int,
    counts: int,
    iterations: int = 100,
    min_fl_count: int = 2,
    step: int = 10 ** 4,
):
    sizes = list(range(0, total + 1, step))
    print(f"min fl count: {min_fl_count}")
    print("size", "min", "max", "mean", "sd")
    for s in sizes:
        tmp = []
        for i in range(iterations):
            tally = defaultdict(lambda: 0)
            for k in random.sample(counts, s):
                tally[k] += 1
            tmp.append(
                len(tally)
            )  # tmp.append(len(filter(lambda k: tally[k]>=min_fl_count, tally)))
        # tmp = [len(set(random.sample(counts, s))) for i in xrange(iter)]
        _mean = sum(tmp) * 1.0 / len(tmp)
        _std = math.sqrt(sum((x - _mean) ** 2 for x in tmp) * 1.0 / len(tmp))
        logger.info(s, min(tmp), max(tmp), _mean, _std)


@app.command(name="")
def main(
    count_filename: str = typer.Argument(...),
    by: str = typer.Option("id", help="Unique specifier name"),
    iterations: int = typer.Option(100, help="Number of iterations (default: 100)"),
    range: Optional[Tuple[int, int]] = typer.Option(
        None, help="Length range (ex: (1000,2000))"
    ),
    min_fl_count: int = typer.Option(1, help="Minimum FL count"),
    step: int = typer.Option(10000, help="Step size (default: 10000)"),
):
    min_len, max_len = None, None
    if range is not None:
        min_len, max_len = eval(range)
        assert 0 <= min_len < max_len

    total, counts = get_counts(count_filename, min_fl_count, by, min_len, max_len)
    subsample(total, counts, iterations, min_fl_count, step)


if __name__ == "__main__":
    typer.run(main)
