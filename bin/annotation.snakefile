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

rule index_genome:
    input:
        "annotation/genome.fa"
    output:
        "annotation/genome.fa.fai"
    shell:
        "samtools faidx {input}"

rule get_annotation:
    output:
        "annotation/annotation.gtf"
    shell:
        " wget -q -O \
            - {ANNOTATION_URL} \
            | gunzip > {output}"

rule make_txdb:
    input:
        "annotation/annotation.gtf"
    output:
        "annotation/annotation_txdb.sqlite"
    script:
        "txdb.R"

rule sort_annotation:
    input:
        "annotation/annotation.gtf"
    output:
        "annotation/annotation_sort.gtf"
    shell:
        "bedtools sort -i {input} > {output}"
