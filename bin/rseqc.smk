rule rseqc_getBed:
# https://www.biostars.org/p/299573/
# from Devon Ryan
    input:
        "data/annotation/annotation.gtf"
    output:
        temp("data/annotation/annotation.bed")
    shell:
        """
        awk '{{if ($$3 != "gene") print $$0;}}' {input} \
            | grep -v '^#' \
            | gtfToGenePred /dev/stdin /dev/stdout \
            | genePredToBed stdin {output}
        """


rule rseqc_genebodyCoverage_illumina:
    input:
        "data/bam/illumina/{sample}_Aligned.sortedByCoord.out.bam.bai",
        bam = "data/bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bed = "data/annotation/annotation.bed"
    output:
        "data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.curves.pdf",
        "data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.r",
        "data/comparisons/geneBody_coverage/illumina/{sample}.geneBodyCoverage.txt"
    params:
        outputDir = "data/comparisons/geneBody_coverage/illumina/{sample}"
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
        temp("{file}_primaryOnly.bam")
    threads: 4
    shell:
        "samtools view -b -F 0x904 -@ {threads} {input} > {output}"


rule rseqc_genebodyCoverage:
    input:
        "data/bam/{dataset}/merged/{dataset}_{barcode}_genome_primaryOnly.bam.bai",
        bam = "data/bam/{dataset}/merged/{dataset}_{barcode}_genome_primaryOnly.bam",
        bed = "data/annotation/annotation.bed"
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
