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
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "res/comparisons/comparisons_feature_detection.html"
    script:
        "../scripts/comparisons_featureDetection.Rmd"
