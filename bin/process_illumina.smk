# fastqc rules
rule fastqc:
    input:
        fastq = "data/fastq/illumina/{type}/{sample}_R1_001.fastq.gz"
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
        fastq_fw = "data/fastq/illumina/raw/{sample}_R1_001.fastq.gz",
        fastq_rv = "data/fastq/illumina/raw/{sample}_R2_001.fastq.gz"
    output:
        fastq_fw = "data/fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "data/fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz",
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


# STAR rules
rule star_index:
    input:
        genome = "data/annotation/genome.fa",
        annotation = "data/annotation/annotation.gtf"
    output: "indices/STAR/Genome"
    params:
        outputDir = "indices/STAR"
    threads: 20
    shell:
        "STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.outputDir} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.annotation}"


rule star_align:
    input:
        genome = "indices/STAR/Genome",
        fastq_fw = "data/fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "data/fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "data/bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        "data/bam/illumina/{sample}_SJ.out.tab"
    params:
        prefix = "data/bam/illumina/{sample}_",
        inputDir = "indices/STAR"
    threads: 20
    shell:
        "STAR \
            --runMode alignReads \
            --genomeDir {params.inputDir} \
            --readFilesIn {input.fastq_fw} {input.fastq_rv} \
            --genomeLoad LoadAndKeep \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 15000000000 \
            --outFilterType BySJout \
            --outFilterMultimapNmax 20 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 1 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000"


# salmon rules
BEDTOOLS = config["BEDTOOLS"]
MASHMAP = config["MASHMAP"]

rule salmon_prepareDecoys:
    input:
        transcripts = "data/annotation/transcripts.fa",
        genome = "data/annotation/genome.fa",
        annotation = "data/annotation/annotation.gtf"
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
        fastq_fw = "data/fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "data/fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "salmon/{sample}/quant.sf"
    params:
        inputDir = "indices/salmon",
        outputDir = "data/salmon/{sample}"
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


rule sam_to_sortBam:
    input:
        "{file}.sam"
    output:
        "{file}_sort.bam"
    threads:
        5
    shell:
        "samtools sort -l 5 -o {output} -O bam -@ {threads} {input}"


rule samtools_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"


rule make_deseqDataSet_illumina:
    input:
        salmon_out = expand("data/quantification/salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        txdb = "data/annotation/annotation_txdb.sqlite",
        sample_info = config["SAMPLE_INFO"]
    params:
        type = "{type}"
    wildcard_constraints:
        type = "gene|transcript"
    output:
        "data/deseq/illumina/illumina_dds_{type}.rds"
    script:
        "deseq_importIllumina.R"
