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
