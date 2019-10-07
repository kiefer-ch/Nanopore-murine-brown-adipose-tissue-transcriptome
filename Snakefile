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

rule render_deseq_all:
    input:
        "res/txlevel_all/deseq_txlevel_all.html",
        "res/txlevel_ont/deseq_txlevel_ont.html",
        "res/genelevel_ont/deseq_genelevel_ont.html"

rule render_correlations:
    input:
        "res/txlevel_ont/txlevel_ont_cm_rld.csv.gz",
        "res/wien/6samples/ChrKiefer_6samples_reg_log_transf_counts.tsv",
        "res/txlevel_ont/txlevel_ont_cm_tpm.csv.gz",
        "res/genelevel_ont/genelevel_ont_cm_rld.csv.gz",
        "res/genelevel_ont/genelevel_ont_cm_tpm.csv.gz",
        "res/wien/6samples/ChrKiefer_6samples_raw_gene_counts.tsv",
        "sample_info/sampleInfo.csv"
    output:
        "res/comparisons/correlations.html"
    script:
        "bin/comparions_ont_illumina.Rmd"

rule render_GOcomp:
    input:
        "res/txlevel_ont/txlevel_ont_de.csv.gz",
        "res/genelevel_ont/genelevel_ont_de.csv.gz",
        "res/wien/6samples/ChrKiefer_6samples_raw_gene_counts.tsv",
        "sample_info/sampleInfo.csv"
    output:
        "res/comparisons/go.html"
    script:
        "bin/comparions_go_illumina.Rmd"

rule render_ont_dtu:
    threads: 4
    input:
        "data/scaledTPM.rds"
    output:
        "res/dtu_ont/ont_dtu"
    script:
        "bin/drimseq_ont.Rmd"

rule render_all_dtu:
    threads: 4
    input:
        "data/scaledTPM.rds"
    output:
        "res/dtu_ont/ont_dtu"
    script:
        "bin/drimseq_ont.Rmd"
