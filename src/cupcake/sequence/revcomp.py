from typing import List

import typer
from Bio import Seq

from cupcake import version_callback

app = typer.Typer(name="cupcake.sequence.revcomp")


@app.command(name="")
def main(
    sequences: List[str] = typer.Argument(..., help="A list of sequences to convert"),
    version: bool = typer.Option(
        None,
        "--version",
        callback=version_callback,
        is_eager=True,
        help="Prints the version of the SQANTI3 package.",
    ),
) -> None:
    [print(Seq.Seq(seq).reverse_complement()) for seq in sequences]


if __name__ == "__main__":
    typer.run(main)
