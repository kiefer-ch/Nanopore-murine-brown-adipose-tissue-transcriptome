import pandas as pd

# URLS to annotation
GENCODE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22"

GENOME_URL = expand(
    "{base_url}/GRCm38.primary_assembly.genome.fa.gz", base_url=GENCODE_URL)
TRANSCRIPTS_URL = expand(
    "{base_url}/gencode.vM22.transcripts.fa.gz", base_url=GENCODE_URL)
ANNOTATION_URL = expand(
    "{base_url}/gencode.vM22.primary_assembly.annotation.gtf.gz", base_url=GENCODE_URL)

# Sample IDs
SAMPLES = pd.read_csv("sample_info/illumina_sample_id.txt", sep='\t')["sample"].tolist()

# packrat rules
rule packrat_init:
    script:
        "bin/packrat_init.R"

# annotation rules
rule get_transcripts:
    output:
        "annotation/transcripts.fa"
    shell:
        " wget -q -O \
            - {TRANSCRIPTS_URL} \
            | gunzip > {output}"

rule get_genome:
    output:
        "annotation/genome.fa"
    shell:
        " wget -q -O \
            - {GENOME_URL} \
            | gunzip > {output}"

rule get_annotation:
    output:
        "annotation/annotation.gtf"
    shell:
        " wget -q -O \
            - {ANNOTATION_URL} \
            | gunzip > {output}"

# STAR rules
rule generate_STARgenome:
    input:
        genome = "annotation/genome.fa",
        annotation = "annotation/annotation.gtf"
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

rule map_STAR:
    input:
        genome = "indices/STAR/Genome",
        fastq_fw = "fastq/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "fastq/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        prefix = "BAM/{sample}_",
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

rule all_STAR:
    input:
        expand("fastq/trimmed/{sample}_R1_001_trimmed.fastq.gz", sample=SAMPLES),
        expand("fastq/trimmed/{sample}_R2_001_trimmed.fastq.gz", sample=SAMPLES),
        expand("BAM/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        "indices/STAR/Genome"
    params:
        genome = "indices/STAR"
    shell:
        "STAR \
            --runMode alignReads \
            --genomeDir {params.genome} \
            --genomeLoad Remove"

# Make bigwigs
rule index_BAM:
    input:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule make_bigwigs:
    input:
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "BAM/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        fw = "BW/{sample}_fw.bw",
        rv = "BW/{sample}_rv.bw"
    threads: 10
    shell:
        "bamCoverage \
            -b {input.bam} \
            -o {output.fw} \
            --filterRNAstrand forward \
            -p {threads} \
            --effectiveGenomeSize 2652783500  \
            --normalizeUsing BPM && \
        bamCoverage \
            -b {input.bam} \
            -o {output.rv} \
            --filterRNAstrand reverse \
            -p {threads} \
            --effectiveGenomeSize 2652783500 \
            --normalizeUsing BPM"

rule make_hub:
    input:
        expand("BW/{sample}_fw.bw", sample=SAMPLES),
        expand("BW/{sample}_rv.bw", sample=SAMPLES),
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "nanoporeibat_hub/hub/hub.txt",
        "nanoporeibat_hub/hub/genomes.txt",
        "nanoporeibat_hub/hub/mm10/trackDb.txt",
        expand("nanoporeibat_hub/BW/{sample}_fw.bw", sample=SAMPLES),
        expand("nanoporeibat_hub/BW/{sample}_rv.bw", sample=SAMPLES)
    params:
       url = "http://bioinformatik.sdu.dk/solexa/webshare/christoph/nanoporeibat_hub"
    script:
        "bin/generateUCSChub.R {input.sample_info} {params.url}"

