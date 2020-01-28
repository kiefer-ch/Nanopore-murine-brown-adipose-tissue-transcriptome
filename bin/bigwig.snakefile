# Make bigwigs
rule samtools_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"

rule bamCoverage_stranded:
    input:
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
        bai = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam.bai"
    output:
        fw = "bw/illumina/{sample}_fw.bw",
        rv = "bw_illumina/{sample}_rv.bw"
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

rule merge_bam_teloprime:
    input:
        bam1 = "bam/teloprime/X1_flowcell/20191107_X1_{type}_{barcode}_q7_sort.bam",
        bai1 = "bam/teloprime/X1_flowcell/20191107_X1_{type}_{barcode}_q7_sort.bam.bai",
        bam2 = "bam/teloprime/X3_flowcell/20191107_X3_{type}_{barcode}_q7_sort.bam",
        bai2 = "bam/teloprime/X3_flowcell/20191107_X3_{type}_{barcode}_q7_sort.bam.bai"
    wildcard_constraints:
        type = "genome|transcriptome"
    output:
        "bam/teloprime/{barcode}_{type}.bam"
    shell:
        "samtools merge {output} {input.bam1} {input.bam2}"


rule merge_bam_cDNA:
    input:
        bam1 = "bam/cdna/pool1/20200108_pool1_{type}_{barcode}_q7_sort.bam",
        bai1 = "bam/cdna/pool1/20200108_pool1_{type}_{barcode}_q7_sort.bam.bai",
        bam2 = "bam/cdna/pool2/20200108_pool2_{type}_{barcode}_q7_sort.bam",
        bai2 = "bam/cdna/pool2/20200108_pool2_{type}_{barcode}_q7_sort.bam.bai"
    wildcard_constraints:
        type = "genome|transcriptome"
    output:
        "bam/cdna/{barcode}_{type}.bam"
    shell:
        "samtools merge {output} {input.bam1} {input.bam2}"


rule bamCoverage_nonstranded:
    input:
        bam = "bam/teloprime/{barcode}_genome.bam",
        bai = "bam/teloprime/{barcode}_genome.bam.bai"
    output:
        "bw/teloprime/{barcode}.bw"
    threads: 10
    shell:
        "bamCoverage \
            -b {input.bam} \
            -o {output} \
            -p {threads} \
            --effectiveGenomeSize 2652783500  \
            --normalizeUsing BPM"

rule makeHub:
    input:
        expand("bw/illumina/{sample}_fw.bw", sample=SAMPLES),
        expand("bw/illumina/{sample}_rv.bw", sample=SAMPLES),
        expand("bw/teloprime/{barcode}.bw", barcode=BARCODES),
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
