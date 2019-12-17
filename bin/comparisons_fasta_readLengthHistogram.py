#!/usr/bin/env python3

"""fasta_tx_lengths.py: get transcript lenghts from a transcriptome fasta."""

__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"

from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import numpy as np

count = 0
lengths = []

with open(snakemake.input[0]) as in_handle:
    for title, seq in SimpleFastaParser(in_handle):
        count += 1
        lengths.append(len(seq))

values, counts = np.unique(lengths, return_counts=True)
data = {'read_length': values, "n": counts}
df = pd.DataFrame.from_records(data)

df.to_csv(snakemake.output[0], index=False, na_rep='0')
