#!/usr/bin/env python3

"""bam_countFlags.py: Count flags from multiple bams into a dataframe."""

__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"

import pysam
import pandas as pd
import numpy as np

df_list = []

for file in snakemake.input:

    print(file)

    unmapped = 0
    primary = 0
    supplementary = 0

    inFile = pysam.AlignmentFile(file)

    for read in inFile:
        if read.flag == 4:
            unmapped += 1
        elif read.flag == 0 or read.flag == 16:
            primary += 1
        elif read.flag == 256 or read.flag == 275:
            supplementary += 1

    data = data = {"unmapped": unmapped, "primary": primary, "supplementary": supplementary, "file": file}
    df_list.append(data)

df = pd.DataFrame(df_list)
df.to_csv(snakemake.output[0], index=False, na_rep='0')
