#!/usr/bin/env python3

"""fastq_readLengthHistogram.py: Aggregate read lengths for multiple fastq files."""

__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip
import pandas as pd
import numpy as np

df = pd.DataFrame({'read_length': []})

for file in snakemake.input:
    print(file)
    count = 0
    lengths = []

    with gzip.open(file, "rt") as in_handle:
        for title, seq, qual in FastqGeneralIterator(in_handle):
            count += 1
            lengths.append(len(seq))

    values, counts = np.unique(lengths, return_counts=True)
    data = {'read_length': values, file: counts}
    df_ = pd.DataFrame.from_records(data)
    df = pd.merge(df, df_, on='read_length', how='outer')

df.to_csv(snakemake.output[0], index=False, na_rep='0')
