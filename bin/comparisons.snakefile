# comparisons
rule quantification_correlation:
    input:
        teloprime_tx_raw = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_cts.csv.gz",
        teloprime_tx_norm = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_ntd.csv.gz",
        illumina_tx_raw = "res/deseq/illumina/txlevel_ont/txlevel_ont_cm_cts.csv.gz",
        tpm_tx = "res/deseq/illumina/txlevel_ont/txlevel_ont_cm_tpm.csv.gz",
        teloprime_gene_raw = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_cts.csv.gz",
        teloprime_gene_norm = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_ntd.csv.gz",
        illumina_gene_raw = "res/deseq/illumina/genelevel_ont/genelevel_ont_cm_cts.csv.gz",
        tpm_gene = "res/deseq/illumina/genelevel_ont/genelevel_ont_cm_tpm.csv.gz",
        sample_info = "sample_info/sampleInfo.csv",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biotype_groups = "data/biotype_groups.csv"
    output:
        "res/comparisons/quantification_correlation.html"
    script:
        "quantification_correlation.Rmd"

rule feature_detection:
    input:
        teloprime_tx = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_cts.csv.gz",
        illumina_tx = "res/deseq/illumina/txlevel_ont/txlevel_ont_cm_cts.csv.gz",
        tpm_tx = "res/deseq/illumina/txlevel_ont/txlevel_ont_cm_tpm.csv.gz",
        teloprime_gene = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_cts.csv.gz",
        illumina_gene = "res/deseq/illumina/genelevel_ont/genelevel_ont_cm_cts.csv.gz",
        sample_info = "sample_info/sampleInfo.csv",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv"
    output:
        "res/comparisons/feature_detection.html"
    script:
        "feature_detection.Rmd"

rule count_reads_bam:
    input:
        "BAM/bam_ont/{file}.bam",
        "BAM/bam_ont/{file}.bam.bai",
    output:
        "data/countReads/{file}.rds"
    script:
        "count_reads.R"

rule count_reads_fastq_teloprime:
    input:
        expand(
            "fastq/teloprime/X{flowcell}_flowcell/{barcode}_q7.fastq.gz", barcode=BARCODES, flowcell=[1, 3])
    output:
        "data/countReads/teloprime_fastqReadLengths.csv"
    script:
        "fastq_counts.py"

rule read_lengths:
    input:
        expand("data/countReads/{barcode}.rds", barcode=BARCODES),
        expand("data/countReads/{barcode}_transcriptome.rds", barcode=BARCODES),
        teloprime_readLengths = "data/countReads/teloprime_fastqReadLengths.csv",
        annotation_txLengths = "data/annotation/transcripts_lengths.csv"
    output:
        "res/comparisons/read_lengths.html"
    script:
        "read_lengths.Rmd"

rule count_txLength_reference:
    input:
        "annotation/transcripts.fa"
    output:
        "data/annotation/transcripts_lengths.csv"
    script:
        "fasta_tx_lengths.py"

rule compare_differentialExpressionAnalysis:
    input:
        illumina_tx = "res/deseq/illumina/txlevel_ont/txlevel_ont_de.csv.gz",
        illumina_gene = "res/deseq/illumina/genelevel_ont/genelevel_ont_de.csv.gz",
        illumina_dtu = "res/drimseq/illumina/dtu_ont/stageR_drimseq_dtu.csv.gz",
        teloprime_tx = "res/deseq/teloprime/txlevel/teloprime_txlevel_de.csv.gz",
        teloprime_gene = "res/deseq/teloprime/genelevel/teloprime_genelevel_de.csv.gz",
        teloprime_dtu = "res/wien/DRIMSeq_stageR/stageR/stageR_final_output_padj_GeneSymbols.tsv",
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
    output:
        "res/comparisons/comparisons_dgeDteDtu_onlyDetectedByBoth.html"
    script:
        "comparisons_dgeDte_onlyDetectedByBoth.Rmd"

rule browser_plots:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        illumina_bam = "BAM/5035_S34_Aligned.sortedByCoord.out.bam",
        illumina_bam_warm = "BAM/5034_S33_Aligned.sortedByCoord.out.bam",
        teloprime_bam = "BAM/bam_ont/barcode02.bam",
        teloprime_bam_warm = "BAM/bam_ont/barcode01.bam"
    params:
        bw_folder = "BW"
    output:
        adcy3 = "res/browser_plots/adcy3.pdf",
        adcy3_bam = "res/browser_plots/adcy3_bam.pdf",
        adcy3_bam_ausschnitt = "res/browser_plots/adcy3_bam_ausschnitt.pdf",
        ctcflos = "res/browser_plots/ctcflos.pdf",
        ctcflos_bam = "res/browser_plots/ctcflos_bam.pdf",
        gm15551 = "res/browser_plots/Gm15551.pdf",
        gm15551_bam = "res/browser_plots/Gm15551_bam.pdf"
    script:
        "browser_plots.R"


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
