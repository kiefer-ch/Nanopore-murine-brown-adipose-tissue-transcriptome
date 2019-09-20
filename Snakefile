__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"

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
    run:
        shell("""
            R -e 'source("packrat/init.R")'
        """)
        shell("R -e 'packrat::restore()'")

# Include other rules
include: "bin/annotation.snakefile"

include: "bin/fastq.snakefile"

include: "bin/star.snakefile"

include: "bin/bigwig.snakefile"

include: "bin/rseqc.snakefile"

include: "bin/salmon.snakefile"

include: "bin/txImport.snakefile"

# render rmd rules
rule render_ont_gene:
    input:
        "data/dds_gencode.vM22_gene_ont.rds"
    output:
        "res/genelevel_ont/deseq_genelevel_ont.html"
    script:
        "bin/deseq_genelevel_ont.Rmd"

rule render_all_gene:
    input:
        "data/dds_gencode.vM22_gene.rds"
    output:
        "res/genelevel_ont/deseq_genelevel_all.html"
    script:
        "bin/deseq_genelevel_all.Rmd"

rule render_ont_transcript:
    input:
        "data/dds_gencode.vM22_transcript_ont.rds"
    output:
        "res/txlevel_ont/deseq_txlevel_ont.html"
    script:
        "bin/deseq_txlevel_ont.Rmd"

rule render_all_transcript:
    input:
        "data/dds_gencode.vM22_transcript.rds"
    output:
        "res/txlevel_all/deseq_txlevel_all.html"
    script:
        "bin/deseq_txlevel_all.Rmd"
