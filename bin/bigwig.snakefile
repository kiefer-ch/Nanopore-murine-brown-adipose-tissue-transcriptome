# Make bigwigs
rule samtools_index:
    input:
        "{file}.bam"
    output:
        "{file}.bam.bai"
    shell:
        "samtools index {input}"


def get_bamnames_bw(wildcards):
    if wildcards.dataset == "cdna":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["cdna"]
        filename = "bam/{}/{}_genome.bam".format(wildcards.dataset, barcode)
    elif wildcards.dataset == "teloprime":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["ont"]
        filename = "bam/{}/{}_genome.bam".format(wildcards.dataset, barcode)
    elif wildcards.dataset == "illumina":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["illumina"]
        filename = "bam/illumina/{}_Aligned.sortedByCoord.out.bam".format(barcode)
    return filename


def get_bainames_bw(wildcards):
    if wildcards.dataset == "cdna":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["cdna"]
        filename = "bam/{}/{}_genome.bam.bai".format(wildcards.dataset, barcode)
    elif wildcards.dataset == "teloprime":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["ont"]
        filename = "bam/{}/{}_genome.bam.bai".format(wildcards.dataset, barcode)
    elif wildcards.dataset == "illumina":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["illumina"]
        filename = "bam/illumina/{}_Aligned.sortedByCoord.out.bam.bai".format(barcode)
    return filename


rule bamCoverage_stranded:
    input:
        bam = get_bamnames_bw,
        bai = get_bainames_bw
    output:
        fw = "bw/{dataset}/{sample}_{dataset}_fw.bw",
        rv = "bw/{dataset}/{sample}_{dataset}_rv.bw"
    wildcard_constraints:
        dataset = "illumina"
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
        bam = get_bamnames_bw,
        bai = get_bainames_bw
    output:
        "bw/{dataset}/{sample}_{dataset}.bw"
    threads: 10
    wildcard_constraints:
        dataset = "cdna|teloprime"
    shell:
        "bamCoverage \
            -b {input.bam} \
            -o {output} \
            -p {threads} \
            --effectiveGenomeSize 2652783500  \
            --normalizeUsing BPM"


rule gtfToGenePred:
    input:
        "flair/{dataset}/flair.collapse.isoforms.gtf"
    output:
        temp("nanoporeibat_hub/bigGenePred/flair_{dataset}.genePred")
    wildcard_constraints:
        dataset = "cdna|teloprime"
    shell:
        "gtfToGenePred -genePredExt {input} {output}"


rule genePredToBigGenePred:
    input:
        "{file}.genePred"
    output:
        temp("{file}.txt")
    shell:
        "genePredToBigGenePred {input} {output}"


rule getBigGenePredHelper:
    output:
        "data/bigGenePred.as"
    shell:
        "wget -q -O \
            - https://genome.ucsc.edu/goldenPath/help/examples/bigGenePred.as \
            > {output}"

rule bedSort:
    input:
        "{file}.txt"
    output:
        temp("{file}_sort.txt")
    shell:
        "bedSort {input} {output}"

rule bed2bb:
    input:
        txt = "{file}_sort.txt",
        chromSizes = "annotation/genome.fa.fai",
        helper = "data/bigGenePred.as"
    output:
        "{file}.bb"
    shell:
        "bedToBigBed -type=bed12+8 -tab \
            -as={input.helper} \
            {input.txt} \
            {input.chromSizes} \
            {output}"


rule makeHub:
    input:
        expand("bw/illumina/{sample}_illumina_fw.bw",
            sample=SAMPLE_INFO_ont.index),
        expand("bw/illumina/{sample}_illumina_rv.bw",
            sample=SAMPLE_INFO_ont.index),
        expand("bw/{dataset}/{sample}_{dataset}.bw",
            dataset=["teloprime", "cdna"], sample=SAMPLE_INFO_ont.index),
        expand("nanoporeibat_hub/bigGenePred/flair_{dataset}.bigGenePred",
            dataset=["teloprime", "cdna"]),
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "nanoporeibat_hub/hub/hub.txt",
        "nanoporeibat_hub/hub/genomes.txt",
        "nanoporeibat_hub/hub/mm10/trackDb.txt",
        expand("nanoporeibat_hub/bw/{sample}_illumina_fw.bw", sample=SAMPLES),
        expand("nanoporeibat_hub/bw/{sample}_illumina_rv.bw", sample=SAMPLES),
        expand("nanoporeibat_hub/bw/{sample}_{dataset}_rv.bw",
            dataset=["teloprime", "cdna"], sample=SAMPLES)
    params:
        url = "http://bioinformatik.sdu.dk/solexa/webshare/christoph/nanoporeibat_hub"
    script:
        "UCSChub.R"
