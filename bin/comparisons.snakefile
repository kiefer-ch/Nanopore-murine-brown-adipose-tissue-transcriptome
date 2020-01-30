# comparisons
rule quantification_correlation:
    input:
        tx_counts = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_ntd.csv.gz",
            dataset=["cdna", "teloprime"]),
        tx_tpm = "res/deseq/illumina/txlevel/illumina_txlevel_cm_ntd.csv.gz",
        gene_counts = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_cm_ntd.csv.gz",
            dataset=["cdna", "teloprime"]),
        gene_tpm = "res/deseq/illumina/genelevel/illumina_genelevel_cm_ntd.csv.gz",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv"
    output:
        "res/comparisons/comparisons_quantification_correlation.html"
    script:
        "comparisons_quantification_correlation.Rmd"


rule feature_detection:
    input:
        tx_counts = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
            dataset=["cdna", "teloprime", "illumina"]),
        gene_counts = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_cm_cts.csv.gz",
            dataset=["cdna", "teloprime", "illumina"]),
        gene_tpm = "res/deseq/illumina/genelevel/illumina_genelevel_cm_ntd.csv.gz",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        fig_folder = "res/fig/feature_detection"
    output:
        "res/comparisons/comparisons_feature_detection.html"
    script:
        "comparisons_featureDetection.Rmd"


rule bam_getAlignedLength:
    input:
        "bam/{dataset}/{file}_{type}.bam",
        "bam/{dataset}/{file}_{type}.bam.bai"
    wildcard_constraints:
        type = "genome|transcriptome",
        datatset = "teloprime"
    output:
        "res/comparisons/countReads/{dataset}_{file}_bam_{type}.rds"
    script:
        "comparisons_bam_getAlignedLength.R"


rule bam_getCoverage:
    input:
        "bam/{dataset}/{file}_transcriptome.bam",
        "bam/{dataset}/{file}_transcriptome.bam.bai"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    output:
        "res/comparisons/coverage/{dataset}/{file}.rds"
    script:
        "comparisons_bam_getCoverage.R"


rule coverage:
    input:
        geneBodyCoverage_teloprime = expand("res/comparisons/geneBody_coverage/teloprime/{barcode}.geneBodyCoverage.txt",
            barcode=SAMPLE_INFO_ont["ont"]),
        geneBodyCoverage_cdna =  expand("res/comparisons/geneBody_coverage/cdna/{barcode}.geneBodyCoverage.txt",
            barcode=SAMPLE_INFO_ont["cdna"]),
        geneBodyCoverage_illumina = expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.txt",
            sample=SAMPLES_ont),
        coverage_teloprime = expand("res/comparisons/coverage/teloprime/{barcode}.rds",
            barcode=SAMPLE_INFO_ont["ont"]),
        coverage_cdna = expand("res/comparisons/coverage/cdna/{barcode}.rds",
            barcode=SAMPLE_INFO_ont["cdna"]),
        sample_info = "sample_info/sampleInfo.csv",
        biomaRt_tx = "annotation/biomaRt_tx.rds"
    params:
        fig_folder = "res/fig/coverage"
    output:
        "res/comparisons/comparisons_coverage.html"
    script:
        "comparisons_coverage.Rmd"


rule fastq_readLengthHistogram:
    input:
        expand(
            "fastq/teloprime/X{flowcell}_flowcell/{barcode}_q7.fastq.gz", barcode=BARCODES, flowcell=[1, 3])
    output:
        "res/comparisons/countReads/teloprime_fastqReadLengths.csv"
    script:
        "comparisons_fastq_readLengthHistogram.py"


rule fasta_readLengthHistogram:
    input:
        "annotation/transcripts.fa"
    output:
        "annotation/annotation_transcript_lengths.csv"
    script:
        "comparisons_fasta_readLengthHistogram.py"


rule read_lengths:
    input:
        teloprime_bam_genome = expand("res/comparisons/countReads/teloprime_{barcode}_bam_genome.rds", barcode=BARCODES),
        teloprime_bam_tx = expand("res/comparisons/countReads/teloprime_{barcode}_bam_transcriptome.rds", barcode=BARCODES),
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        teloprime_readLengths = "res/comparisons/countReads/teloprime_fastqReadLengths.csv",
        annotation_txLengths = "annotation/annotation_transcript_lengths.csv"
    params:
        fig_folder = "res/fig/read_lengths"
    output:
        "res/comparisons/comparisons_readLengths.html"
    script:
        "comparisons_read_lengths.Rmd"


rule compare_differentialExpressionAnalysis:
    input:
        illumina_tx = "res/deseq/illumina/txlevel_ont/txlevel_ont_de.csv.gz",
        illumina_gene = "res/deseq/illumina/genelevel_ont/genelevel_ont_de.csv.gz",
        illumina_dtu = "res/drimseq/illumina/dtu_ont/stageR_drimseq_dtu.csv.gz",
        teloprime_tx = "res/deseq/teloprime/txlevel/teloprime_txlevel_de.csv.gz",
        teloprime_gene = "res/deseq/teloprime/genelevel/teloprime_genelevel_de.csv.gz",
        teloprime_dtu = "res/wien/teloprime_old/DRIMSeq_stageR/stageR/stageR_final_output_padj_GeneSymbols.tsv",
        teloprime_counts = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_cts.csv.gz",
        illumina_counts = "res/deseq/illumina/genelevel_ont/genelevel_ont_cm_cts.csv.gz",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        grouped_biotypes = "data/biotype_groups.csv",
        tpm_gene = "res/deseq/illumina/genelevel_ont/genelevel_ont_cm_tpm.csv.gz"
    output:
        "res/comparisons/comparisons_dgeDteDtu.html"
    script:
        "comparisons_dgeDteDtu.Rmd"


rule compare_differentialExpressionAnalysis2:
    input:
        illumina_tx = "res/deseq/illumina/txlevel_ont/illumina_txlevel_ont_dds.rds",
        illumina_gene = "res/deseq/illumina/genelevel_ont/illumina_genelevel_ont_dds.rds",
        teloprime_tx = "res/deseq/teloprime/txlevel/teloprime_txlevel_dds.rds",
        teloprime_gene = "res/deseq/teloprime/genelevel/teloprime_genelevel_dds.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        grouped_biotypes = "data/biotype_groups.csv"
    threads: 4
    output:
        "res/comparisons/comparisons_dgeDte_onlyDetectedByBoth.html"
    script:
        "comparisons_dgeDte_onlyDetectedByBoth.Rmd"


rule compare_dtu:
    input:
        illumina_drim = "res/drimseq/illumina/illumina_drimSeqStageR.csv",
        teloprime_drim = "res/drimseq/teloprime/teloprime_drimSeqStageR.csv",
        teloprime_drim_flair = "res/drimseq/teloprime_flair/teloprime_flair_drimSeqStageR.csv",
        illumina_dex = "res/dexseq/illumina/illumina_dexseq_results.csv.gz",
        teloprime_dex = "res/dexseq/teloprime/teloprime_dexseq_results.csv.gz",
        chip = "res/chip/gencode_full_range.csv.gz",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    output:
        "res/comparisons/comparisons_dtu.html"
    script:
        "comparisons_dtu.Rmd"





rule render_GOcomp:
    input:
        "res/txlevel_ont/txlevel_ont_de.csv.gz",
        "res/genelevel_ont/genelevel_ont_de.csv.gz",
        "res/wien/6samples/ChrKiefer_6samples_raw_gene_counts.tsv",
        "sample_info/sampleInfo.csv"
    output:
        "res/comparisons/go.html"
    script:
        "bin/comparions_go_ont.Rmd"
