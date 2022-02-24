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
    output:
        "res/comparisons/comparisons_coverage.html"
    script:
        "comparisons_coverage.Rmd"


# unnormalised quantification correlation
rule quantification_getAverageCountsTables:
    input:
        ont_quant = [expand("data/quantification/teloprime/merged/teloprime_merged_{barcode}_quant.tsv", barcode=SAMPLE_INFO_ont["ont"]),
                     expand("data/quantification/cdna/merged/cdna_merged_{barcode}_quant.tsv", barcode=SAMPLE_INFO_ont["cdna"]),
                     expand("data/quantification/rna/merged/rna_merged_{barcode}_quant.tsv", barcode=["rt", "cool"])],
        illumina_quant = expand("data/quantification/salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        txdb = "data/annotation/annotation_txdb.sqlite"
    output:
        gene = "data/comparisons/counts/comparisons_meanCounts_gene.tsv.gz",
        tx = "data/comparisons/counts/comparisons_meanCounts_tx.tsv.gz"
    script:
        "comparisons_quantification_averageCountsTable.R"


rule quantification_correlation:
    input:
        biomaRt_gene = "data/annotation/biomaRt_gene.rds",
        biomaRt_tx = "data/annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        gene = "data/comparisons/counts/comparisons_meanCounts_gene.tsv.gz",
        tx = "data/comparisons/counts/comparisons_meanCounts_tx.tsv.gz"
    output:
        "res/comparisons/comparisons_quantification_correlation.html"
    script:
        "comparisons_quantification_correlation.Rmd"


rule quantification_getCountsTables:
    input:
        ont_quant = [expand("data/quantification/teloprime/{flowcell}/teloprime_{flowcell}_{barcode}_quant.tsv",
                         flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["ont"]),
                     expand("data/quantification/cdna/{flowcell}/cdna_{flowcell}_{barcode}_quant.tsv",
                         flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["cdna"])],
    output:
        "data/comparisons/counts/comparisons_counts.tsv.gz"
    script:
        "comparisons_quantification_countsTable.R"


rule quantification_correlationWithinSamples:
    input:
        "data/comparisons/counts/comparisons_counts.tsv.gz"
    output:
        "res/comparisons/comparisons_quantification_correlationWithinSamples.html"
    script:
        "comparisons_quantification_correlationWithinSamples.Rmd"


# normalised quantification correlation
rule quantification_correlation_normalised:
    input:
        biomaRt_gene = "data/annotation/biomaRt_gene.rds",
        biomaRt_tx = "data/annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        gene = [expand("data/deseq/{dataset}/{dataset}_gene_cm_ntd.csv.gz",
                    dataset=["cdna", "rna", "teloprime"]),
                 "data/deseq/illumina/illumina_gene_cm_tpm.csv.gz"],
        tx = [expand("data/deseq/{dataset}/{dataset}_transcript_cm_ntd.csv.gz",
                  dataset=["cdna", "rna", "teloprime"]),
              "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz"]
    output:
        "res/comparisons/comparisons_quantification_correlation_normalised.html"
    script:
        "comparisons_quantification_correlation_normalised.Rmd"


rule make_deseqDataSet_libraries:
    input:
        teloprime = expand("data/quantification/teloprime/{flowcell}/teloprime_{flowcell}_{barcode}_quant.tsv",
            flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["ont"]),
        cdna = expand("data/quantification/cdna/{flowcell}/cdna_{flowcell}_{barcode}_quant.tsv",
            flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["cdna"]),
        sample_info = config["SAMPLE_INFO"]
    output:
        "data/comparisons/counts/counts_withinSamplesNorm.tsv.gz"
    script:
        "comparisons_quantification_importONT_library.R"


rule quantification_correlationWithinSamples_normalised:
    input:
        "data/comparisons/counts/counts_withinSamplesNorm.tsv.gz"
    output:
        "res/comparisons/comparisons_quantification_correlationWithinSamples_normalised.html"
    script:
        "comparisons_quantification_correlationWithinSamples_normalised.Rmd"


rule counts_pca:
    input:
        gene = [expand("data/deseq/{dataset}/{dataset}_gene_cm_ntd.csv.gz",
                    dataset=["cdna", "rna", "teloprime"]),
                 "data/deseq/illumina/illumina_gene_cm_tpm.csv.gz"],
        tx = [expand("data/deseq/{dataset}/{dataset}_transcript_cm_ntd.csv.gz",
                  dataset=["cdna", "rna", "teloprime"]),
              "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz"],
        sample_info = config["SAMPLE_INFO"]
    output:
        "res/comparisons/comparisons_countsPCA.html"
    script:
        "comparisons_countsPCA.Rmd"


# feature_detection
rule feature_detection:
    input:
        tx_counts = expand("data/deseq/{dataset}/{dataset}_transcript_cm_cts.csv.gz",
            dataset=["cdna", "teloprime", "illumina", "rna"]),
        gene_counts = expand("data/deseq/{dataset}/{dataset}_gene_cm_cts.csv.gz",
            dataset=["cdna", "teloprime", "illumina", "rna"]),
        gene_tpm = "data/deseq/illumina/illumina_gene_cm_tpm.csv.gz",
        tx_tpm = "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz",
        biomaRt_gene = "data/annotation/biomaRt_gene.rds",
        biomaRt_tx = "data/annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        sample_info = config["SAMPLE_INFO"]
    params:
        cutoff = 1
    output:
        "res/comparisons/comparisons_feature_detection.html"
    script:
        "comparisons_featureDetection.Rmd"


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
    script:
        "comparisons_dgeDte.Rmd"


# reannotation
rule compare_reannotation:
    input:
        gffcompare =
            [expand("data/reannotation/flair/annotation/gffcmp.{dataset}_flair.isoforms.gtf.tmap",
                dataset=["cdna", "teloprime", "rna"]),
            expand("data/reannotation/stringtie/gffcmp.{dataset}_stringtie.gtf.tmap",
                dataset=["cdna", "teloprime", "rna", "illumina"])],
        sqanti =
            [expand("data/comparisons/reannotation/squanti/stringtie/{dataset}/{dataset}_stringtie_noUnknownStrand_classification.txt",
                dataset=["cdna", "teloprime", "rna", "illumina"]),
            expand("data/comparisons/reannotation/squanti/flair/{dataset}/{dataset}_flair.isoforms_classification.txt",
                dataset=["cdna", "teloprime", "rna"])],
        gffcompare_stats = "data/comparisons/reannotation/gffcompare/RQ/gffcmp.stats"
    output:
        "res/comparisons/comparisons_reannotation.html"
    script:
        "comparisons_reannotation.Rmd"


rule compare_dtu:
    input:
        dtu_res = expand("data/drimseq/{dataset}_{file}_dtu_res.csv",
            dataset = ["illumina", "cdna"],
            file = ["ref", "illumina_stringtie", "cdna_flair", "teloprime_stringtie"]),
        qpcr = "data/qpcr/dtu_qpcr.csv"
    output:
        "res/comparisons/comparisons_dtu.html",
    params:
        dtu_cutoff = 0.1
    script:
        "comparisons_dtu.Rmd"
