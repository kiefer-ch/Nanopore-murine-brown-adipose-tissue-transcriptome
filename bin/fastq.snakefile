# fastqc rules
rule run_fastqc_raw:
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

rule qc_raw_all:
    input:
        expand("qc/fastqc/raw/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/raw/{sample}_R2_001_fastqc.zip", sample=SAMPLES)
    output:
        "qc/multiqc_raw.html",
        "qc/multiqc_raw_data.zip"
    shell:
        "multiqc -f -z \
            qc/fastqc/raw \
            -o qc \
            -n multiqc_raw.html"


rule run_fastqc_trimmed:
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

rule qc_trimmed_all:
    input:
        expand("qc/fastqc/raw/{sample}_R1_001_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/raw/{sample}_R2_001_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/trimmed/{sample}_R1_001_trimmed_fastqc.zip", sample=SAMPLES),
        expand("qc/fastqc/trimmed/{sample}_R2_001_trimmed_fastqc.zip", sample=SAMPLES)
    output:
        "qc/multiqc_trimmed.html"
        "qc/multiqc_trimmed_data.zip"
    shell:
        "multiqc -f -z \
            -c bin/.multiqc.conf \
            qc/fastqc \
            -o qc \
            -n multiqc_trimmed.html"

# cutadapt rules
rule trim_adapters:
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
