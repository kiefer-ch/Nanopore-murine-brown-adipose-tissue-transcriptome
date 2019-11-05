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

rule count_reads:
    input:
        "BAM/bam_ont/{file}.bam",
        "BAM/bam_ont/{file}.bam.bai",
    output:
        "data/countReads/{file}.rds"
    script:
        "count_reads.R"

rule read_lengths:
    input:
        expand("data/countReads/{barcode}.rds", barcode=BARCODES),
        expand("data/countReads/{barcode}_transcriptome.rds", barcode=BARCODES)

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
