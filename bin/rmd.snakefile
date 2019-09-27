# render rmd rules
rule render_ont_gene:
    input:
        "data/dds_gencode.vM22_gene_ont.rds"
    output:
        "res/genelevel_ont/deseq_genelevel_ont.html"
    script:
        "bin/deseq_genelevel_ont.Rmd"

rule render_all_gene:
    input:
        "data/dds_gencode.vM22_gene.rds"
    output:
        "res/genelevel_ont/deseq_genelevel_all.html"
    script:
        "bin/deseq_genelevel_all.Rmd"

rule render_ont_transcript:
    input:
        "data/dds_gencode.vM22_transcript_ont.rds"
    output:
        "res/txlevel_ont/deseq_txlevel_ont.html"
    script:
        "bin/deseq_txlevel_ont.Rmd"

rule render_all_transcript:
    input:
        "data/dds_gencode.vM22_transcript.rds"
    output:
        "res/txlevel_all/deseq_txlevel_all.html"
    script:
        "bin/deseq_txlevel_all.Rmd"

rule render_deseq_all:
    input:
        "res/txlevel_all/deseq_txlevel_all.html",
        "res/txlevel_ont/deseq_txlevel_ont.html",
        "res/genelevel_ont/deseq_genelevel_ont.html"

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
        "bin/comparions_go_illumina.Rmd"

