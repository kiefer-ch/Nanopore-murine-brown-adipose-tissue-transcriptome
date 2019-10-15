# fastqc rules
rule fastqc_raw:
    input:
        fastq_fw = "fastq/raw/{sample}_R1_001.fastq.gz",
        fastq_rv = "fastq/raw/{sample}_R2_001.fastq.gz"
    output:
        "qc/fastqc/raw/{sample}_R1_001_fastqc.html",
        "qc/fastqc/raw/{sample}_R1_001_fastqc.zip",
        "qc/fastqc/raw/{sample}_R2_001_fastqc.html",
        "qc/fastqc/raw/{sample}_R2_001_fastqc.zip"
    params:
        outputDir = "qc/fastqc/raw"
    shell:
        "fastqc {input.fastq_fw} \
            --noextract \
            --quiet \
            -o {params.outputDir} && \
        fastqc {input.fastq_rv} \
            --noextract \
            --quiet \
            -o {params.outputDir}"

# cutadapt rules
rule cutadapt_trim:
    input:
        fastq_fw = "fastq/raw/{sample}_R1_001.fastq.gz",
        fastq_rv = "fastq/raw/{sample}_R2_001.fastq.gz"
    output:
        fastq_fw = "fastq/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "fastq/trimmed/{sample}_R2_001_trimmed.fastq.gz"
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
            {input.fastq_fw} {input.fastq_rv}"

rule fastqc_trimmed:
    input:
        fastq_fw = "fastq/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "fastq/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "qc/fastqc/trimmed/{sample}_R1_001_trimmed_fastqc.html",
        "qc/fastqc/trimmed/{sample}_R1_001_trimmed_fastqc.zip",
        "qc/fastqc/trimmed/{sample}_R2_001_trimmed_fastqc.html",
        "qc/fastqc/trimmed/{sample}_R2_001_trimmed_fastqc.zip"
    params:
        outputDir = "qc/fastqc/trimmed"
    shell:
        "fastqc {input.fastq_fw} \
            --noextract \
            --quiet \
            -o {params.outputDir} && \
        fastqc {input.fastq_rv} \
            --noextract \
            --quiet \
            -o {params.outputDir}"
