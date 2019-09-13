# salmon rules
rule prepare_decoys:
    input:
        transcripts = "annotation/transcripts.fa",
        genome = "annotation/genome.fa",
        annotation = "annotation/annotation.gtf"
    output:
        "indices/salmon/decoy/gentrome.fa",
        "indices/salmon/decoy/decoys.txt"
    params:
        outputDir = "indices/salmon/decoy"
    threads: 20
    shell:
        "bin/generateDecoyTranscriptome.sh \
            -j {threads} \
            -g {input.genome} \
            -t {input.transcripts} \
            -a {input.annotation} \
            -b {BEDTOOLS} \
            -m {MASHMAP} \
            -o {params.outputDir}"

rule generate_salmonIndex:
    input:
        gentrome = "indices/salmon/decoy/gentrome.fa",
        decoy = "indices/salmon/decoy/decoys.txt"
    output:
        "indices/salmon/duplicate_clusters.tsv",
        "indices/salmon/hash.bin",
        "indices/salmon/rsd.bin",
        "indices/salmon/sa.bin",
        "indices/salmon/txpInfo.bin"
    params:
        outputDir = "indices/salmon"
    threads: 20
    shell:
        "salmon index \
            --gencode \
            -t {input.gentrome} \
            -i {params.outputDir} \
            -d {input.decoy} \
            -p {threads}"

rule map_salmon:
    threads: 5
    input:
        "indices/salmon/duplicate_clusters.tsv",
        "indices/salmon/hash.bin",
        "indices/salmon/rsd.bin",
        "indices/salmon/sa.bin",
        "indices/salmon/txpInfo.bin",
        fastq_fw = "fastq/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "fastq/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "salmon/{sample}/quant.sf"
    params:
        inputDir = "indices/salmon",
        outputDir = "salmon/{sample}"
    shell:
        "salmon quant \
            --gcBias \
            --seqBias \
            --numGibbsSamples 25 \
            -i {params.inputDir} \
            -l A \
            -1 {input.fastq_fw} \
            -2 {input.fastq_rv} \
            -p {threads} \
            --validateMappings \
            -o {params.outputDir}"

rule all_salmon:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES)

rule qc_salmon:
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
        expand("qc/RSeQC/junction_saturation/{sample}.junctionSaturation_plot.r", sample=SAMPLES),
        expand("salmon/{sample}/quant.sf", sample=SAMPLES)
    output:
        "qc/multiqc_salmon.html",
        "qc/multiqc_salmon_data.zip"
    shell:
        "multiqc -f -z \
            -c bin/.multiqc.conf \
            qc/fastqc/ BAM/ salmon/ qc/RSeQC \
            -o qc \
            -n multiqc_salmon"
