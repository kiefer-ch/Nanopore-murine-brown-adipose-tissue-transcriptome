# fastqc rules
rule fastqc:
    input:
        fastq = "fastq/illumina/{type}/{sample}_R1_001.fastq.gz"
    output:
        "qc/fastqc/{type}/{sample}_fastqc.html",
        "qc/fastqc/{type}/{sample}_fastqc.zip"
    params:
        outputDir = "qc/fastqc/{type}"
    shell:
        "fastqc {input.fastq} \
            --noextract \
            --quiet \
            -o {params.outputDir}"

# cutadapt rules
rule cutadapt_trim:
    input:
        fastq_fw = "fastq/illumina/raw/{sample}_R1_001.fastq.gz",
        fastq_rv = "fastq/illumina/raw/{sample}_R2_001.fastq.gz"
    output:
        fastq_fw = "fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz",
        report = "qc/cutadapt/{sample}.txt"
    threads: 5
    shell:
        "cutadapt \
            -j {threads} \
            -q 28 \
            -m 30 \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            -o {output.fastq_fw} \
            -p {output.fastq_rv} \
            {input.fastq_fw} {input.fastq_rv} \
            > {output.report}"
