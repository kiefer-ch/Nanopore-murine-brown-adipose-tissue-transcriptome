# DGE
rule deseq_genelevel_ont:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/genelevel_ont"
    input:
        "data/dds_gencode.vM22_gene_ont.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/genelevel_ont/deseq_genelevel_ont.html"
    script:
        "deseq_genelevel_ont.Rmd"

rule deseq_genelevel_all:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/genelevel_all"
    input:
        "data/dds_gencode.vM22_gene.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/genelevel_all/deseq_genelevel_all.html"
    script:
        "deseq_genelevel_all.Rmd"

# DTE
rule deseq_txlevel_ont:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/txlevel_ont"
    input:
        "data/dds_gencode.vM22_transcript_ont.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/txlevel_ont/deseq_txlevel_ont.html"
    script:
        "deseq_txlevel_ont.Rmd"

rule deseq_txlevel_all:
    threads: 4
    params:
        out_folder = "res/deseq/illumina/txlevel_all"
    input:
        "data/dds_gencode.vM22_transcript.rds",
        "data/biotype_groups.csv"
    output:
        "res/deseq/illumina/txlevel_all/deseq_txlevel_all.html"
    script:
        "deseq_txlevel_all.Rmd"

# teloprime
rule deseq_teloprime:
    input:
        ont_gene_raw = "res/wien/6samples/ChrKiefer_6samples_raw_gene_counts.tsv",
        sample_info = "sample_info/sampleInfo.csv",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    threads: 4
    output:
        gene_rld = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_rld.csv.gz",
        gene_cts = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_cts.csv.gz",
        gene_ntd = "res/deseq/teloprime/genelevel/teloprime_genelevel_cm_ntd.csv.gz",
        gene_de = "res/deseq/teloprime/genelevel/teloprime_genelevel_de.csv.gz",
        tx_ntd = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_ntd.csv.gz",
        tx_rld = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_rld.csv.gz",
        tx_cts = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_cts.csv.gz",
        tx_de = "res/deseq/teloprime/txlevel/teloprime_txlevel_de.csv.gz",
    script:
        "deseq_teloprime.R"
