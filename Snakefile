__author__ = "Christoph Kiefer"
__email__ = "christophak@bmb.sdu.dk"


import pandas as pd


configfile: "config.yaml"


# URLS to annotation
GENCODE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22"

GENOME_URL = expand(
    "{base_url}/GRCm38.primary_assembly.genome.fa.gz", base_url=GENCODE_URL)
TRANSCRIPTS_URL = expand(
    "{base_url}/gencode.vM22.transcripts.fa.gz", base_url=GENCODE_URL)
ANNOTATION_URL = expand(
    "{base_url}/gencode.vM22.primary_assembly.annotation.gtf.gz", base_url=GENCODE_URL)


# Sample IDs
SAMPLE_INFO = pd.read_csv("sample_info/sampleInfo.csv", sep=',')
SAMPLE_INFO = SAMPLE_INFO.set_index("sample_id")
SAMPLE_INFO_illumina = SAMPLE_INFO[SAMPLE_INFO["illumina"].notnull()]
SAMPLES = SAMPLE_INFO_illumina["illumina"].tolist()
SAMPLE_INFO_ont = SAMPLE_INFO[SAMPLE_INFO["ont"].notnull()]
SAMPLES_ont = SAMPLE_INFO_ont["illumina"].tolist()
BARCODES = SAMPLE_INFO_illumina[SAMPLE_INFO_illumina["ont"].notnull()]["ont"].tolist()


# target rules
rule all:
    input:
        expand("BW/bw_ont/{barcode}.bw", barcode=BARCODES),
        expand("BW/{sample}_fw.bw", sample=SAMPLES),
        expand("BW/{sample}_rv.bw", sample=SAMPLES),
        "qc/multiqc_aligned.html",
        "res/deseq/illumina/txlevel_ont/deseq_txlevel_ont.html",
        "res/deseq/illumina/genelevel_ont/deseq_genelevel_ont.html",
        "res/deseq/illumina/txlevel_all/deseq_txlevel_all.html",
        "res/deseq/illumina/genelevel_all/deseq_genelevel_all.html",
        "res/dtu/dtu_all/all_dtu.html",
        "res/dtu/dtu_ont/ont_dtu.html",
        "res/comparisons/feature_detection.html",
        "res/comparisons/quantification_correlation.html"


rule illumina_align:
    input:
        expand("qc/RSeQC/bam_stat/{sample}", sample=SAMPLES),
        expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf", sample=SAMPLES),
        expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.r", sample=SAMPLES),
        expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.DupRate_plot.pdf", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.DupRate_plot.r", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.pos.DupRate.xls", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.seq.DupRate.xls", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.junction.bed", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.junction.xls", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.junction_plot.r", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.splice_events.pdf", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.splice_junction.pdf", sample=SAMPLES),
        expand("qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.pdf", sample=SAMPLES),
        expand("qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.r", sample=SAMPLES),
        expand("salmon/{sample}/quant.sf", sample=SAMPLES)
    output:
        "qc/multiqc_aligned.html",
        "qc/multiqc_aligned_data.zip"
    shell:
        "multiqc -f -z \
            -c bin/.multiqc.conf \
            qc/fastqc/ qc/cutadapt bam/illumina/ salmon/ qc/RSeQC \
            -o qc \
            -n multiqc_aligned"


rule illumina_trimm:
    input:
        expand("fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz", sample=SAMPLES),
        expand("fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz", sample=SAMPLES),
        expand("qc/fastqc/raw/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/raw/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/trimmed/{sample}_R1_001_trimmed_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/trimmed/{sample}_R2_001_trimmed_fastqc.zip", sample=SAMPLES),
    output:
        "qc/multiqc_trimmed.html",
        "qc/multiqc_trimmed_data.zip"
    shell:
        "multiqc -f -z \
            -c bin/multiqc.conf \
            qc/fastqc qc/cutadapt \
            -o qc \
            -n multiqc_trimmed.html"

rule dexseq_all:
    input:
        expand("res/dexseq/{dataset}/{dataset}_heatmap.html",
               dataset=["illumina", "teloprime", "cdna"]),
        expand("res/dexseq/{dataset}/{dataset}_dexseq_results.csv.gz",
               dataset=["illumina", "teloprime", "cdna"])


rule drimseq_all:
    input:
        expand("res/drimseq/{dataset}/{dataset}_drimSeqStageR.html",
               dataset=["illumina", "teloprime", "cdna", "cdna_flair", "teloprime_flair"])


rule comparisons_all:
    input:
        "res/comparisons/comparisons_quantification_correlation.html",
        "res/comparisons/comparisons_countsPCA.html",
        "res/comparisons/comparisons_feature_detection.html",
        "res/comparisons/comparisons_coverage.html",
        "res/comparisons/comparisons_dtu.html",
        "res/comparisons/comparisons_readLengths_fastq.html",
        "res/comparisons/comparisons_readLengths_bam.html",
        "res/comparisons/comparisons_dgeDte.html",
        "res/comparisons/comparisons_reannotation.html"


rule flair_all:
    input :
        expand("flair/{dataset}/flair.collapse.{dataset}.isoforms.gtf",
            dataset=["teloprime", "cdna", "rna"])


# Include other rules
include: "bin/annotation.snakefile"
include: "bin/fastq.snakefile"
include: "bin/star.snakefile"
include: "bin/bigwig.snakefile"
include: "bin/rseqc.snakefile"
include: "bin/salmon.snakefile"
include: "bin/comparisons.smk"
include: "bin/dexseq.snakefile"
include: "bin/deseq.snakefile"
include: "bin/drimseq.snakefile"
include: "bin/flair.snakefile"
include: "bin/stringtie.snakefile"
include: "bin/gffcompare.snakefile"
include: "bin/qpcr.snakefile"
include: "bin/browserTracks.smk"
