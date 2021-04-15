from typing import List

import typer
from Bio import Seq

app = typer.Typer(name="cupcake.sequence.revcomp")


@app.command(name="")
def main(
    sequences: List[str] = typer.Argument(..., help="A list of sequences to convert")
) -> None:
    [print(Seq.Seq(seq).reverse_complement()) for seq in sequences]


if __name__ == "__main__":
    typer.run(main)
