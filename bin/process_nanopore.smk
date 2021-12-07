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





rule minimap_mapGenome:
    input:
        index = "indices/minimap2/genome_minimap2.mmi",
        fastq = get_minimapInput
    output:
        sam = "data/bam/{dataset}/{dataset}_{flowcell}_genome_{barcode}_q7.sam"
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
        index = "indices/minimap2/transcriptome_minimap2.mmi",
        fastq = get_minimapInput
    output:
        sam = "data/bam/{dataset}/{dataset}_{flowcell}_transcriptome_{barcode}_q7.sam"
    threads: 24
    shell:
        "minimap2 \
            -ax map-ont \
            -t {threads} \
            --secondary=no \
            -uf \
            {input.index} \
            {input.fastq} > {output.sam}"


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


def get_deseqONTinput(wildcards):
    if wildcards.dataset == "teloprime":
        files = expand("data/quantification/teloprime/merged/teloprime_merged_{barcode}_quant.tsv", barcode=SAMPLE_INFO_ont["ont"])
    elif wildcards.dataset == "cdna":
        files = expand("data/quantification/cdna/merged/cdna_merged_{barcode}_quant.tsv", barcode=SAMPLE_INFO_ont["cdna"])
    elif wildcards.dataset == "rna":
        files = expand("data/quantification/rna/merged/rna_merged_{barcode}_quant.tsv", barcode=["cool", "rt"])
    return files


rule make_deseqDataSet_ONT:
    input:
        counts = get_deseqONTinput,
        txdb = "data/annotation/annotation_txdb.sqlite",
        sample_info = config["SAMPLE_INFO"]
    params:
        type = "{type}"
    wildcard_constraints:
        type = "gene|transcript",
        dataset = "cdna|rna|teloprime"
    output:
        "data/deseq/{dataset}/{dataset}_dds_{type}.rds"
    script:
        "deseq_importONT.R"
