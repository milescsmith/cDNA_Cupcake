#!/usr/bin/env python
import math
import random
from collections import defaultdict
from csv import DictReader
from typing import List, Tuple

import typer

from cupcake.logging import cupcake_logger as logger

app = typer.Typer(name="cupcake.annotation.subsample_with_category")


def get_counts(
    count_filename: str, min_fl_count=2, key="id", min_len=None, max_len=None
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
            count_d[f"{r[key]}---{r['category']}"] += c
            total += c

    counts = []
    for k, v in count_d.items():
        counts += [k] * v

    return total, counts


def subsample(
    total: int, counts: List[int], iter=100, min_fl_count=2, step=10 ** 4
) -> None:
    sizes = list(range(100, total + 1, step))
    logger.info(f"min fl count: {min_fl_count}")
    logger.info("size\tcategory\tmin\tmax\tmean\tsd")
    for s in sizes:
        tmp = defaultdict(lambda: [])  # category --> N iterations
        for i in range(iter):
            tally = defaultdict(lambda: 0)
            uniq_id_count = defaultdict(lambda: set())  # category -> unique ids
            for k in random.sample(counts, s):
                tally[k] += 1
            for k in tally:
                _id, _cat = k.split("---")
                uniq_id_count[_cat].add(
                    _id
                )  # if tally[k] >= min_fl_count: uniq_id_count[_cat].add(_id)
            for _cat in uniq_id_count:
                tmp[_cat].append(len(uniq_id_count[_cat]))

        for _cat in tmp:
            _mean = sum(tmp[_cat]) * 1.0 / len(tmp[_cat])
            _std = math.sqrt(
                sum((x - _mean) ** 2 for x in tmp[_cat]) * 1.0 / len(tmp[_cat])
            )
            logger.info(
                f"{s}\t{_cat}\t{min(tmp[_cat])}\t{max(tmp[_cat])}\t{_mean}\t{_std}"
            )


@app.command(name="")
def main(
    count_filename: str = typer.Argument(...),
    by: str = typer.Option("id", help="Unique specifier name"),
    iterations: int = typer.Option(100, help="Number of iterations"),
    len_range: int = typer.Option(
        None, "--range", help="Length range (ex: (1000,2000), default None)"
    ),
    min_fl_count: int = typer.Option(1, help="Minimum FL count"),
    step: int = typer.Option(10000, help="Step size"),
) -> None:

    min_len, max_len = None, None
    if len_range is not None:
        min_len, max_len = eval(len_range)
        assert 0 <= min_len < max_len

    total, counts = get_counts(count_filename, min_fl_count, by, min_len, max_len)
    subsample(total, counts, iterations, min_fl_count, step)


if __name__ == "__main__":
    typer.run(main)
