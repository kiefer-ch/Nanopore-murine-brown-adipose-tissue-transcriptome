# DGE
rule deseq_genelevel_ont:
    threads: 4
    params:
        out_folder = "res/genelevel_ont"
    input:
        "data/dds_gencode.vM22_gene_ont.rds",
        "data/biotype_groups.csv"
    output:
        "res/genelevel_ont/deseq_genelevel_ont.html"
    script:
        "deseq_genelevel_ont.Rmd"

rule deseq_genelevel_all:
    threads: 4
    params:
        out_folder = "res/genelevel_all"
    input:
        "data/dds_gencode.vM22_gene.rds",
        "data/biotype_groups.csv"
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
        "data/dds_gencode.vM22_transcript_ont.rds",
        "data/biotype_groups.csv"
    output:
        "res/txlevel_ont/deseq_txlevel_ont.html"
    script:
        "deseq_txlevel_ont.Rmd"

rule deseq_txlevel_all:
    threads: 4
    params:
        out_folder = "res/txlevel_all"
    input:
        "data/dds_gencode.vM22_transcript.rds",
        "data/biotype_groups.csv"
    output:
        "res/txlevel_all/deseq_txlevel_all.html"
    script:
        "deseq_txlevel_all.Rmd"
