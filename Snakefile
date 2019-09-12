import pandas as pd

# URLS to annotation
GENCODE_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22"

GENOME_URL = expand(
    "{base_url}/GRCm38.primary_assembly.genome.fa.gz", base_url=GENCODE_URL)
TRANSCRIPTS_URL = expand(
    "{base_url}/gencode.vM22.transcripts.fa.gz", base_url=GENCODE_URL)
ANNOTATION_URL = expand(
    "{base_url}/gencode.vM22.primary_assembly.annotation.gtf.gz", base_url=GENCODE_URL)

# PAths to software required by salmon
BEDTOOLS = "/data/bin/bedtools/bedtools-v2.28/bedtools"
MASHMAP = "/data/home/christophak/bin/mashmap"

# Sample IDs
SAMPLES = pd.read_csv("sample_info/illumina_sample_id.txt", sep='\t', header=None)[0].tolist()

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

rule annotation_all:
    input:
        "annotation/transcripts.fa"
        "annotation/genome.fa"
        "annotation/annotation.gtf"

# fastq rules
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
        expand("BAM/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)
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
        expand("nanoporeibat_hub/bw/{sample}_fw.bw", sample=SAMPLES),
        expand("nanoporeibat_hub/bw/{sample}_rv.bw", sample=SAMPLES)
    params:
        url = "http://bioinformatik.sdu.dk/solexa/webshare/christoph/nanoporeibat_hub"
    shell:
        "bin/generateUCSChub.R {input.sample_info} {params.url}"

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

rule tximport_gene:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_gene.rds"
    shell:
        "bin/txImport_gene.R {input.sample_info} {input.annotation} {output} FALSE"

rule tximport_transcript:
    input:
        expand("salmon/{sample}/quant.sf", sample=SAMPLES),
        annotation = "annotation/annotation.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dds_gencode.vM22_transcript.rds"
    shell:
        "bin/txImport_gene.R {input.sample_info} {input.annotation} {output} TRUE"

rule tximport_all:
    input:
        "data/dds_gencode.vM22_transcript.rds",
        "data/dds_gencode.vM22_gene.rds"
