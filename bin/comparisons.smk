# fastq read characterisation
def get_fastqnames(wildcards):
    files = list()
    if wildcards.dataset == "teloprime":
        for barcode in SAMPLE_INFO_ont["ont"]:
            for flowcell in ["X1_flowcell", "X3_flowcell"]:
                filename = "fastq/{}/{}/{}_q7.fastq.gz".format(
                    wildcards.dataset, flowcell, barcode)
                files.append(filename)
    elif wildcards.dataset == "cdna":
        for barcode in SAMPLE_INFO_ont["cdna"]:
            for pool in ["pool1", "pool2"]:
                filename = "fastq/{}/{}/{}_q7.fastq.gz".format(
                    wildcards.dataset, pool, barcode)
                files.append(filename)
    elif wildcards.dataset == "rna":
        for temperature in ["rt", "cool"]:
            filename = "fastq/{}/{}_q7.fastq.gz".format(
                wildcards.dataset, temperature)
            files.append(filename)
    return files


rule readLength_fastq_histogram:
    input:
        get_fastqnames
    output:
        "data/comparisons/fastq/readLengthDistribution/{dataset}_fastqReadLengths.csv"
    script:
        "comparisons_read_lengths_fastq_readLengthHistogram.py"


rule readQuality_fastq_histogram:
    input:
        get_fastqnames
    output:
        "data/comparisons/fastq/readQualityDistribution/{dataset}_fastqQualities.csv"
    script:
        "comparisons_read_lengths_fastq_readLengthHistogram.py"


rule readLengths_fastq:
    input:
        readLengths = expand("data/comparisons/fastq/readLengthDistribution/{dataset}_fastqReadLengths.csv",
                             dataset=["teloprime", "cdna", "rna"]),
        readQualities = expand("data/comparisons/fastq/readQualityDistribution/{dataset}_fastqQualities.csv",
                             dataset=["teloprime", "cdna", "rna"]),
        sample_info = config["SAMPLE_INFO"]
    output:
        "res/comparisons/comparisons_readLengths_fastq.html"
    script:
        "comparisons_read_lengths_fastq.Rmd"


# Mapping characterisation
def get_input_bam_AlignedLength(wildcards):
    if wildcards.dataset == "teloprime":
        file_name = "data/bam/teloprime/merged/teloprime_{}_{}.bam".format(
            wildcards.barcode, wildcards.type)
    elif wildcards.dataset == "cdna":
        file_name = "data/bam/cdna/merged/cdna_{}_{}.bam".format(
            wildcards.barcode, wildcards.type)
    elif wildcards.dataset == "rna":
        file_name = "data/bam/rna/merged/rna_{}_{}.bam".format(
            wildcards.barcode, wildcards.type)
    return [file_name, file_name + ".bai"]


rule readLengths_bam:
    input:
        get_input_bam_AlignedLength
    wildcard_constraints:
        type = "genome|transcriptome",
        dataset = "teloprime|cdna|rna"
    output:
        "data/comparisons/mapping/{dataset}/{dataset}_{barcode}_bam_{type}.rds"
    script:
        "comparisons_read_lengths_bam.R"


rule readLengths_bam_collapseTranscripts:
    input:
        expand("data/comparisons/mapping/teloprime/teloprime_{barcode}_bam_{{type}}.rds",
            barcode=SAMPLE_INFO_ont["ont"]),
        expand("data/comparisons/mapping/cdna/cdna_{barcode}_bam_{{type}}.rds",
            barcode=SAMPLE_INFO_ont["cdna"]),
        expand("data/comparisons/mapping/rna/rna_{barcode}_bam_{{type}}.rds",
            barcode=["rt", "cool"])
    wildcard_constraints:
        type = "genome|transcriptome"
    output:
        "data/comparisons/mapping/bam_countReads_collapsed_{type}.rds"
    script:
        "comparisons_read_lengths_bam_collapseTranscripts.R"


rule readLengths_bam_report_transcriptome:
    input:
        collapsed_transcripts = "data/comparisons/mapping/bam_countReads_collapsed_transcriptome.rds",
        biomaRt_tx = "data/annotation/biomaRt_tx.rds",
        sample_info = config["SAMPLE_INFO"]
    output:
        "res/comparisons/comparisons_readLengths_bam_transcriptome.html"
    script:
        "comparisons_read_lengths_bam_transcriptome.Rmd"


rule readLengths_bam_report_genome:
    input:
        collapsed_transcripts = "data/comparisons/mapping/bam_countReads_collapsed_genome.rds",
        sample_info = config["SAMPLE_INFO"]
    output:
        "res/comparisons/comparisons_readLengths_bam_genome.html"
    script:
        "comparisons_read_lengths_bam_genome.Rmd"


# Coverage
rule coverage_getCoverage:
    input:
        "data/bam/{dataset}/merged/{dataset}_{barcode}_transcriptome.bam"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    output:
        "data/comparisons/coverage/{dataset}/merged/{dataset}_{barcode}_coverage.rds"
    script:
        "comparisons_coverage_getCoverage.R"


rule coverage_collapse:
    input:
        coverage = [expand("data/comparisons/coverage/teloprime/merged/teloprime_{barcode}_coverage.rds", barcode=SAMPLE_INFO_ont["ont"]),
                    expand("data/comparisons/coverage/cdna/merged/cdna_{barcode}_coverage.rds", barcode=SAMPLE_INFO_ont["cdna"]),
                    expand("data/comparisons/coverage/rna/merged/rna_{barcode}_coverage.rds", barcode=["rt", "cool"])],
        biomaRt_tx = "data/annotation/biomaRt_tx.rds"
    output:
        "data/comparisons/coverage/collapsed_coverage.rds"
    script:
        "comparisons_coverage_collapseCoverage.R"


rule coverage_report:
    input:
        geneBodyCoverage_teloprime = expand("data/comparisons/geneBody_coverage/teloprime/{barcode}.geneBodyCoverage.txt",
            barcode=SAMPLE_INFO_ont["ont"]),
        geneBodyCoverage_cdna = expand("data/comparisons/geneBody_coverage/cdna/{barcode}.geneBodyCoverage.txt",
            barcode=SAMPLE_INFO_ont["cdna"]),
        geneBodyCoverage_illumina = expand("data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.txt",
            sample=SAMPLES_ont),
        geneBodyCoverage_rna = expand("data/comparisons/geneBody_coverage/rna/{barcode}.geneBodyCoverage.txt",
            barcode=["rt", "cool"]),
        collapsed_coverage = "data/comparisons/coverage/collapsed_coverage.rds",
        tpm = "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz",
        sample_info = config["SAMPLE_INFO"],
        biomaRt_tx = "data/annotation/biomaRt_tx.rds"
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "res/comparisons/comparisons_coverage.html"
    script:
        "../scripts/comparisons_coverage.Rmd"


# dge
rule compare_dge:
    input:
        dte = expand("data/deseq/{dataset}/{dataset}_deseqResults_transcript.csv.gz",
                     dataset=["illumina", "teloprime", "cdna"]),
        dge = expand("data/deseq/{dataset}/{dataset}_deseqResults_gene.csv.gz",
                     dataset=["illumina", "teloprime", "cdna"]),
        gene_counts = [expand("data/deseq/{dataset}/{dataset}_gene_cm_ntd.csv.gz",
                    dataset=["cdna", "teloprime"]),
                 "data/deseq/illumina/illumina_gene_cm_tpm.csv.gz"],
        tx_counts = [expand("data/deseq/{dataset}/{dataset}_transcript_cm_ntd.csv.gz",
                  dataset=["cdna", "teloprime"]),
              "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz"],
        grouped_biotypes = "data/biotype_groups.csv"
    output:
        "res/comparisons/comparisons_dgeDte.html"
    params:
        cutoff = .05
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "../scripts/comparisons_dgeDte.Rmd"
