# RSeQC rules
rule run_bamstat:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "qc/RSeQC/bam_stat/{sample}"
    params:
        outputDir = "qc/RSeQC/bam_stat/{sample}"
    shell:
        "bam_stat.py \
            -i {input.bam} \
            > {params.outputDir}"

rule extract_junctionBED:
    input:
        "data/mm10_Gencode_vM20_gencode.bed.gz"
    output:
        "data/mm10_Gencode_vM20_gencode.bed"
    shell:
        "gunzip {input}"

rule run_genebody_coverage:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
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

rule run_read_duplication:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
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

rule run_junction_annotation:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
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

rule run_junction_saturation:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
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

rule qc_bam_all:
    input:
        expand("qc/RSeQC/bam_stat/{sample}", sample=SAMPLES),
        expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf", sample=SAMPLES),
        expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.r", sample=SAMPLES),
        expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.txt", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.DupRate_plot.pdf", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.DupRate_plot.r", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.pos.DupRate.xls", sample=SAMPLES),
        expand("qc/RSeQC/read_duplication/{sample}.seq.DupRate.xls", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.junction.bed", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.junction.xls", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.junction_plot.r", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.splice_events.pdf", sample=SAMPLES),
        expand("qc/RSeQC/junction_annotation/{sample}.splice_junction.pdf", sample=SAMPLES),
        expand("qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.pdf", sample=SAMPLES),
        expand("qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.r", sample=SAMPLES)
    output:
        "qc/multiqc_bam.html",
        "qc/multiqc_bam_data.zip"
    shell:
        "multiqc -f -z \
            -c bin/.multiqc.conf \
            qc/fastqc/ BAM/ qc/RSeQC \
            -o qc \
            -n multiqc_bam"
