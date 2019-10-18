# comparisons
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
        "bin/comparions_go_ont.Rmd"
