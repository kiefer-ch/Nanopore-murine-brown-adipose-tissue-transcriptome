# DGE
rule deseq_genelevel_ont:
    threads: 4
    params:
        out_folder = "res/genelevel_ont"
    input:
        "data/dds_gencode.vM22_gene_ont.rds"
    output:
        "res/genelevel_ont/deseq_genelevel_ont.html"
    script:
        "deseq_genelevel_ont.Rmd"

rule deseq_genelevel_all:
    threads: 4
    params:
        out_folder = "res/genelevel_all"
    input:
        "data/dds_gencode.vM22_gene.rds"
    output:
        "res/genelevel_all/deseq_genelevel_all.html"
    script:
        "deseq_genelevel_all.Rmd"

# DTE
rule deseq_txlevel_ont:
    threads: 4
    params:
        out_folder = "res/txlevel_ont"
    input:
        "data/dds_gencode.vM22_transcript_ont.rds"
    output:
        "res/txlevel_ont/deseq_txlevel_ont.html"
    script:
        "deseq_txlevel_ont.Rmd"

rule deseq_txlevel_all:
    threads: 4
    params:
        out_folder = "res/txlevel_all"
    input:
        "data/dds_gencode.vM22_transcript.rds"
    output:
        "res/txlevel_all/deseq_txlevel_all.html"
    script:
        "deseq_txlevel_all.Rmd"

# DTU
rule render_ont_dtu:
    threads: 4
    input:
        "data/scaledTPM_ont.rds"
    output:
        "res/dtu_ont/ont_dtu.html"
    script:
        "bin/drimseq_ont.Rmd"

rule render_all_dtu:
    threads: 4
    input:
        "data/scaledTPM_all.rds"
    output:
        "res/dtu_all/all_dtu.html"
    script:
        "bin/drimseq_all.Rmd"

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
