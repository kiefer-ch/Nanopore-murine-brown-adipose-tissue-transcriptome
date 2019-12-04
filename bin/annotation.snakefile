rule get_transcripts:
    output:
        "annotation/transcripts.fa"
    shell:
        " wget -q -O \
            - {TRANSCRIPTS_URL} \
            | gunzip > {output}"

rule get_genome:
    output:
        "annotation/genome.fa"
    shell:
        " wget -q -O \
            - {GENOME_URL} \
            | gunzip > {output}"

rule index_fasta:
    input:
        "{genome}.fa"
    output:
        "{genome}.fa.fai"
    shell:
        "samtools faidx {input}"

rule get_annotation:
    output:
        "annotation/annotation.gtf"
    shell:
        " wget -q -O \
            - {ANNOTATION_URL} \
            | gunzip > {output}"

rule make_txdb_from_gtf:
    input:
        "{annotation}.gtf"
    output:
        "{annotation}_txdb.sqlite"
    script:
        "txdb.R"

rule sort_gtf:
    input:
        "{annotation}.gtf"
    output:
        "{annotation}_sort.gtf"
    shell:
        "bedtools sort -i {input} > {output}"

rule get_biomart_tx:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
    output:
        tx = "annotation/biomaRt_tx.rds"
    script:
        "biomaRt_tx.R"

rule get_biomart_gene:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
    output:
        gene = "annotation/biomaRt_gene.rds"
    script:
        "biomaRt_gene.R"
