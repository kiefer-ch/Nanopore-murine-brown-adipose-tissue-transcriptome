rule minimap_index:
    input:
        "data/annotation/{type}.fa"
    output:
        "indices/minimap2/{type}_minimap2.mmi"
    threads: 3
    shell:
        "minimap2 \
        -t {threads} \
        -d {output} \
        {input}"


def get_minimapInput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/fastq/{}/{}/{}_q7.fastq.gz".format(
            wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/fastq/rna/{}_q7.fastq.gz".format(wildcards.barcode)
    return file_name


def get_minimapGenomeOutput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/bam/{}/genome/{}/{}_genome.bam".format(
            wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/bam/rna/genome/{}_genome.bam".format(wildcards.barcode)
    return file_name


def get_minimapTranscriptomeOutput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/bam/{}/transcriptome/{}/{}_genomtranscriptome.bam".format(
            wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/bam/rna/transcriptome/{}_transcriptome.bam".format(wildcards.barcode)
    return file_name


rule minimap_mapGenome:
    input:
        index = "indices/minimap2/genome_minimap2.mmi",
        fastq = get_minimapInput
    output:
        sam = get_minimapGenomeOutput
    threads: 24
    shell:
        "minimap2 \
            -ax splice \
            -t {threads} \
            --secondary=no \
            -uf \
            {input.index} \
            {input.fastq} > {output.sam}"


rule minimap_mapTranscriptome:
    input:
        index = "indices/minimap2/transcriptome_minimap2.mmi"
        fastq = get_minimapInput
    output:
        sam = get_minimapTranscriptomeOutput
    threads: 24
    shell:
        "minimap2 \
            -ax map-ont \
            -t {threads} \
            --secondary=no \
            -uf \
            {input.index} \
            {input.fastq} > {output.sam}"
