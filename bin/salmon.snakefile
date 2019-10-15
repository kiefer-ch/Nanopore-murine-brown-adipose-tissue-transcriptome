# salmon rules
rule salmon_prepareDecoys:
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

rule salmon_index:
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

rule salmon_align:
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
