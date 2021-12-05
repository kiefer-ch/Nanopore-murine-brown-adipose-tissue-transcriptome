# URLS to annotation
GENCODE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_{}".format(
    config["GENCODE_VERSION"])
GENOME_URL = "{}/GRCm38.primary_assembly.genome.fa.gz".format(GENCODE_URL)
TRANSCRIPTS_URL = "{}/gencode.v{}.transcripts.fa.gz".format(
    GENCODE_URL, config["GENCODE_VERSION"])
ANNOTATION_URL = "{}/gencode.v{}.primary_assembly.annotation.gtf.gz".format(
    GENCODE_URL, config["GENCODE_VERSION"])


rule get_transcripts:
    output:
        "data/annotation/transcripts.fa"
    shell:
        " wget -q -O \
            - {TRANSCRIPTS_URL} \
            | gunzip > {output}"


rule get_genome:
    output:
        "data/annotation/genome.fa"
    shell:
        " wget -q -O \
            - {GENOME_URL} \
            | gunzip > {output}"


rule index_fasta:
    input:
        "{file}.fa"
    output:
        "{file}.fa.fai"
    shell:
        "samtools faidx {input}"


rule get_annotation:
    output:
        "data/annotation/annotation.gtf"
    shell:
        " wget -q -O \
            - {ANNOTATION_URL} \
            | gunzip > {output}"


rule make_txdb_from_gtf:
    input:
        "{file}.gtf"
    output:
        "{file}_txdb.sqlite"
    script:
        "annotation_txdb.R"


rule sort_gtf:
    input:
        "{file}.gtf"
    output:
        "{annofiletation}_sort.gtf"
    shell:
        "bedtools sort -i {input} > {output}"


rule get_biomart_tx:
    input:
        txdb = "data/annotation/annotation_txdb.sqlite",
    output:
        tx = "data/annotation/biomaRt_tx.rds"
    script:
        "annotation_biomaRt_tx.R"


rule get_biomart_gene:
    input:
        txdb = "data/annotation/annotation_txdb.sqlite",
    output:
        gene = "data/annotation/biomaRt_gene.rds"
    script:
        "annotation_biomaRt_gene.R"
