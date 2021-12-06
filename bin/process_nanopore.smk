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
        file_name = "data/bam/{}/{}/{}_{}_genome_{}_q7.sam".format(
            wildcards.dataset, wildcards.flowcell, wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/bam/rna/rna_genome_{}_q7.sam".format(wildcards.barcode)
    return file_name


def get_minimapTranscriptomeOutput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/bam/{}/{}/{}_{}_transcriptome_{}_q7.sam".format(
            wildcards.dataset, wildcards.flowcell, wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/bam/rna/rna_transcriptome_{}_q7.sam".format(wildcards.barcode)
    return file_name


#rule minimap_mapGenome:
#    input:
#        index = "indices/minimap2/genome_minimap2.mmi",
#        fastq = get_minimapInput
#    output:
#        sam = Only input files can be specified as functions
#    threads: 24
#    shell:
#        "minimap2 \
#            -ax splice \
#            -t {threads} \
#            --secondary=no \
#            -uf \
#            {input.index} \
#            {input.fastq} > {output.sam}"


#rule minimap_mapTranscriptome:
#    input:
#        index = "indices/minimap2/transcriptome_minimap2.mmi",
#        fastq = get_minimapInput
#    output:
#        sam = get_minimapTranscriptomeOutput
#    threads: 24
#    shell:
#        "minimap2 \
#            -ax map-ont \
#            -t {threads} \
#            --secondary=no \
#            -uf \
#            {input.index} \
#            {input.fastq} > {output.sam}"


rule samtools_merge:
    input:
        "data/bam/{dataset}/flowcell1/{dataset}_flowcell1_transcriptome_{barcode}_q7_sort.bam",
        "data/bam/{dataset}/flowcell1/{dataset}_flowcell1_transcriptome_{barcode}_q7_sort.bam"
    output:
        "data/bam/{dataset}/merged/{dataset}_{barcode}_{type}.bam"
    threads: 4
    wildcard_constraints:
        dataset = "teloprime|cdna",
        type = "genome|transcriptome"
    shell:
        "samtools merge -@ {threads} {output} {input}"


rule move_rnaMerge:
    input:
        "data/bam/rna/rna_{type}_{barcode}_q7_sort.bam"
    output:
        "data/bam/rna/merged/rna_{barcode}_{type}.bam"
    wildcard_constraints:
        type = "genome|transcriptome"


def get_minimapQuantifyInput(wildcards):
    if wildcards.flowcell in ["flowcell1", "flowcell2"]:
        file_name = "data/bam/{}/{}/{}_{}_transcriptome_{}_q7_sort.bam".format(
            wildcards.dataset, wildcards.flowcell, wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.flowcell == "merged":
        file_name = "data/bam/{}/merged/{}_{}_transcriptome.bam".format(wildcards.dataset, wildcards.dataset, wildcards.barcode)
    return file_name


rule quantify_minimap:
    input:
        get_minimapQuantifyInput
    output:
        "data/quantification/{dataset}/{flowcell}/{dataset}_{flowcell}_{barcode}_quant.tsv"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna",
        type = "genome|transcriptome"
    script:
        "quantify_minimap.R"
