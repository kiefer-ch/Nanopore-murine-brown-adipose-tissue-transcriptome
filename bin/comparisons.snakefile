# comparisons
rule correlations:
    input:
        "res/txlevel_ont/deseq_txlevel_ont.html",
        ont_tx_rlog = "res/wien/6samples/ChrKiefer_6samples_reg_log_transf_counts.tsv",
        ont_gene_raw = "res/wien/6samples/ChrKiefer_6samples_raw_gene_counts.tsv",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        out_folder = "res/comparisons/correlations",
        illumina_tx_tpm = "res/txlevel_ont/txlevel_ont_cm_tpm.csv.gz",
        illumina_tx_rlog = "res/txlevel_ont/txlevel_ont_cm_rld.csv.gz",
        illumina_gene_tpm = "res/genelevel_ont/genelevel_ont_cm_tpm.csv.gz",
        illumina_gene_rlog = "res/genelevel_ont/genelevel_ont_cm_rld.csv.gz"
    output:
        "res/comparisons/correlations/correlations.html"
    script:
        "correlations.Rmd"

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
