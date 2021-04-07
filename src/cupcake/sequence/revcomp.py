import sys

from Bio import Seq

for seq in sys.argv[1:]:
    print(Seq.Seq(seq).reverse_complement())
