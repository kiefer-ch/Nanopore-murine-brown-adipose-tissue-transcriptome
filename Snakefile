__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"


import pandas as pd


configfile: "conf/config.yaml"


# Sample IDs
SAMPLE_INFO = pd.read_csv(config["SAMPLE_INFO"], sep=',')
SAMPLE_INFO = SAMPLE_INFO.set_index("sample_id")
SAMPLE_INFO_illumina = SAMPLE_INFO[SAMPLE_INFO["illumina"].notnull()]
SAMPLES = SAMPLE_INFO_illumina["illumina"].tolist()
SAMPLE_INFO_ont = SAMPLE_INFO[SAMPLE_INFO["ont"].notnull()]
SAMPLES_ont = SAMPLE_INFO_ont["illumina"].tolist()
BARCODES = SAMPLE_INFO_illumina[SAMPLE_INFO_illumina["ont"].notnull()]["ont"].tolist()


# Include other rules
include: "bin/annotation.smk"
include: "bin/process_illumina.smk"
#include: "bin/process_nanopore.smk"
include: "bin/rseqc.smk"
include: "bin/comparisons.smk"
include: "bin/dexseq.snakefile"
include: "bin/deseq.snakefile"
include: "bin/drimseq.snakefile"
include: "bin/flair.snakefile"
include: "bin/stringtie.snakefile"
include: "bin/gffcompare.snakefile"
include: "bin/qpcr.snakefile"
include: "bin/browserTracks.smk"
