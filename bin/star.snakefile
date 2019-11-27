# STAR rules
rule star_index:
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

rule star_align:
    input:
        genome = "indices/STAR/Genome",
        fastq_fw = "fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        prefix = "bam/illumina/{sample}_",
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
