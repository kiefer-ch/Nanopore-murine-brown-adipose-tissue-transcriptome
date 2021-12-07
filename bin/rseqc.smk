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
        bam = "data/bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
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


rule samtools_filterPrimary:
    # remove 0x4 = "unmapped", 0x100 = "not primary alingmend" == "secondary",
    # 0x800 = "supplementary"
    input:
        "{file}.bam"
    output:
        temp("{file}_genome_primaryOnly.bam")
    threads: 4
    shell:
        "samtools view -b -F 0x904 -@ {threads} {input} > {output}"


rule rseqc_genebodyCoverage:
    input:
        bam = "data/bam/{dataset}/merged/{dataset}_{barcode}_genome_primaryOnly.bam",
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
