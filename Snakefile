import pandas as pd

# URLS to annotation
GENCODE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22"

GENOME_URL = expand(
    "{base_url}/GRCm38.primary_assembly.genome.fa.gz", base_url=GENCODE_URL)
TRANSCRIPTS_URL = expand(
    "{base_url}/gencode.vM22.transcripts.fa.gz", base_url=GENCODE_URL)
ANNOTATION_URL = expand(
    "{base_url}/gencode.vM22.primary_assembly.annotation.gtf.gz", base_url=GENCODE_URL)

# Paths to software required by salmon
BEDTOOLS = "/data/bin/bedtools/bedtools-v2.28/bedtools"
MASHMAP = "/data/home/christophak/bin/mashmap"

# Sample IDs
SAMPLE_INFO = pd.read_csv("sample_info/sampleInfo.csv", sep=',')
SAMPLE_INFO_illumina = SAMPLE_INFO[SAMPLE_INFO["illumina"].notnull()]
SAMPLES = SAMPLE_INFO_illumina["illumina"].tolist()
SAMPLE_INFO_ont = SAMPLE_INFO[SAMPLE_INFO["ont"].notnull()]
SAMPLES_ont = SAMPLE_INFO_ont["illumina"].tolist()

# packrat rule
rule packrat_init:
    script:
        "bin/packrat_init.R"

# Include other rules
include: "bin/annotation.snakefile"

include: "bin/fastq.snakefile"

include: "bin/star.snakefile"

include: "bin/bigwig.snakefile"

include: "bin/rseqc.snakefile"

include: "bin/salmon.snakefile"

include: "bin/txImport.snakefile"
