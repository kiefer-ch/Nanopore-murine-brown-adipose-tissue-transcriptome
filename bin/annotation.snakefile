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

rule sort_annotation:
    input:
        "annotation/annotation.gtf"
    output:
        "annotation/annotation_sort.gtf"
    shell:
        "bedtools sort -i {input} > {output}"

rule annotation_all:
    input:
        "annotation/transcripts.fa",
        "annotation/genome.fa",
        "annotation/annotation.gtf",
        "annotation/genome.fa.fai",
        "annotation/annotation_sort.gtf"

