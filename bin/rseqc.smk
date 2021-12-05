# RSeQC rules
rule rseqc_junctionBED:
    input:
        "data/mm10_Gencode_vM20_gencode.bed.gz"
    output:
        "data/mm10_Gencode_vM20_gencode.bed"
    shell:
        "gunzip {input}"
#??? where is this file from?

rule rseqc_genebodyCoverage_illumina:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.curves.pdf",
        "data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.r",
        "data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.txt"
    params:
        outputDir = "data/comparisons/geneBody_coverage/{sample}"
    shell:
        "geneBody_coverage.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.outputDir}"


def get_genebodyCoverageInput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/bam/{}/{}_genome.bam".format(
            wildcards.dataset, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/bam/rna/genome_{}_q7_sort.bam".format(wildcards.barcode)
    return file_name


rule rseqc_genebodyCoverage:
    input:
        bam = get_genebodyCoverageInput,
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "data/comparisons/geneBody_coverage/{dataset}/{barcode}.geneBodyCoverage.curves.pdf",
        "data/comparisons/geneBody_coverage/{dataset}/{barcode}.geneBodyCoverage.r",
        "data/comparisons/geneBody_coverage/{dataset}/{barcode}.geneBodyCoverage.txt"
    params:
        outputDir = "data/comparisons/geneBody_coverage/{dataset}/{barcode}"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    shell:
        "geneBody_coverage.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.outputDir}"
