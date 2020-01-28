# RSeQC rules
rule rseqc_bamstat:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "qc/RSeQC/bam_stat/{sample}"
    params:
        outputDir = "qc/RSeQC/bam_stat/{sample}"
    shell:
        "bam_stat.py \
            -i {input.bam} \
            > {params.outputDir}"

rule rseqc_junctionBED:
    input:
        "data/mm10_Gencode_vM20_gencode.bed.gz"
    output:
        "data/mm10_Gencode_vM20_gencode.bed"
    shell:
        "gunzip {input}"


rule rseqc_genebodyCoverage:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf",
        "qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.r",
        "qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.txt"
    params:
        outputDir = "qc/RSeQC/geneBody_coverage/{sample}"
    shell:
        "geneBody_coverage.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.outputDir}"


rule rseqc_genebodyCoverage_ont:
    input:
        bam = "bam/{dataset}/{barcode}_genome.bam",
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "res/comparisons/geneBody_coverage/{dataset}/{barcode}.geneBodyCoverage.curves.pdf",
        "res/comparisons/geneBody_coverage/{dataset}/{barcode}.geneBodyCoverage.r",
        "res/comparisons/geneBody_coverage/{dataset}/{barcode}.geneBodyCoverage.txt"
    params:
        outputDir = "res/comparisons/geneBody_coverage/{dataset}/{barcode}"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    shell:
        "geneBody_coverage.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.outputDir}"


rule rseqc_readDuplication:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "qc/RSeQC/read_duplication/{sample}.DupRate_plot.pdf",
        "qc/RSeQC/read_duplication/{sample}.DupRate_plot.r",
        "qc/RSeQC/read_duplication/{sample}.pos.DupRate.xls",
        "qc/RSeQC/read_duplication/{sample}.seq.DupRate.xls"
    params:
        outputDir = "qc/RSeQC/read_duplication/{sample}"
    shell:
        "read_duplication.py \
            -i {input.bam} \
            -o {params.outputDir}"

rule rseqc_junctionAnnotation:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "qc/RSeQC/junction_annotation/{sample}",
        "qc/RSeQC/junction_annotation/{sample}.junction.bed",
        "qc/RSeQC/junction_annotation/{sample}.junction.xls",
        "qc/RSeQC/junction_annotation/{sample}.junction_plot.r",
        "qc/RSeQC/junction_annotation/{sample}.splice_events.pdf",
        "qc/RSeQC/junction_annotation/{sample}.splice_junction.pdf"
    params:
        outputDir = "qc/RSeQC/junction_annotation/{sample}"
    shell:
        "junction_annotation.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.outputDir} \
            2> {params.outputDir}"

rule rseqc_junctionSaturation:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bed = "data/mm10_Gencode_vM20_gencode.bed"
    output:
        "qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.pdf",
        "qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.r"
    params:
        outputDir = "qc/RSeQC/junction_saturation/{sample}"
    shell:
        "junction_saturation.py \
            -r {input.bed} \
            -i {input.bam} \
            -o {params.outputDir}"
