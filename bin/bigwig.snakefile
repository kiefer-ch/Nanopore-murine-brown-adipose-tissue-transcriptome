# Make bigwigs
rule samtools_index_illumina:
    input:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "BAM/{sample}_Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"

rule bamCoverage_stranded:
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

rule merge_bam_ont_genome:
    input:
        bam1 = "BAM/bam_ont/X1_flowcell/20190808_X1_genome_{barcode}_q7_sort.bam",
        bai1 = "BAM/bam_ont/X1_flowcell/20190808_X1_genome_{barcode}_q7_sort.bam.bai",
        bam2 = "BAM/bam_ont/X3_flowcell/20190808_X3_genome_{barcode}_q7_sort.bam",
        bai2 = "BAM/bam_ont/X3_flowcell/20190808_X3_genome_{barcode}_q7_sort.bam.bai"
    output:
        "BAM/bam_ont/{barcode}.bam"
    shell:
        "samtools merge {output} {input.bam1} {input.bam2}"

rule merge_bam_ont_transcriptome:
    input:
        bam1 = "BAM/bam_ont/X1_flowcell/20190808_X1_transcriptome_{barcode}_q7_sort.bam",
        bai1 = "BAM/bam_ont/X1_flowcell/20190808_X1_transcriptome_{barcode}_q7_sort.bam.bai",
        bam2 = "BAM/bam_ont/X3_flowcell/20190808_X3_transcriptome_{barcode}_q7_sort.bam",
        bai2 = "BAM/bam_ont/X3_flowcell/20190808_X3_transcriptome_{barcode}_q7_sort.bam.bai"
    output:
        "BAM/bam_ont/{barcode}_transcriptome.bam"
    shell:
        "samtools merge {output} {input.bam1} {input.bam2}"

rule samtools_index_ont:
    input:
        "BAM/bam_ont/{barcode}.bam"
    output:
        "BAM/bam_ont/{barcode}.bam.bai"
    shell:
        "samtools index {input}"

rule bamCoverage_nonstranded:
    input:
        bam = "BAM/bam_ont/{barcode}.bam",
        bai = "BAM/bam_ont/{barcode}.bam.bai"
    output:
        "BW/bw_ont/{barcode}.bw"
    threads: 10
    shell:
        "bamCoverage \
            -b {input.bam} \
            -o {output} \
            -p {threads} \
            --effectiveGenomeSize 2652783500  \
            --normalizeUsing BPM"

rule makeHub_illumina:
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
