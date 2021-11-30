#!/usr/bin/env python3

"""fastq_readQualityHistogram.py: Aggregate quality for multiple fastq files."""

__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"

from Bio import SeqIO
import gzip
import pandas as pd
import numpy as np

df = pd.DataFrame({'read_quality': []})

for file in snakemake.input:
    print(file)
    quals = []

    with gzip.open(file, "rt") as in_handle:
        for rec in SeqIO.parse(in_handle, "fastq"):
            quals.append(np.mean(rec.letter_annotations["phred_quality"]))

    values, counts = np.unique(np.round(quals, 2), return_counts=True)
    data = {'read_quality': values, file: counts}
    df_ = pd.DataFrame.from_records(data)
    df = pd.merge(df, df_, on='read_quality', how='outer')

df.to_csv(snakemake.output[0], index=False, na_rep='0')
