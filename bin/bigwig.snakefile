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

rule makeHub:
    input:
        expand("bw/illumina/{sample}_illumina_fw.bw",
            sample=SAMPLE_INFO_ont.index),
        expand("bw/illumina/{sample}_illumina_rv.bw",
            sample=SAMPLE_INFO_ont.index),
        expand("bw/teloprime/{sample}_teloprime.bw",
            sample=SAMPLE_INFO_ont.index),
        expand("bw/cdna/{sample}_cdna.bw",
            sample=SAMPLE_INFO_ont.index),
        sample_info = "sample_info/sampleInfo.csv"
    # output:
    #     "nanoporeibat_hub/hub/hub.txt",
    #     "nanoporeibat_hub/hub/genomes.txt",
    #     "nanoporeibat_hub/hub/mm10/trackDb.txt",
    #     expand("nanoporeibat_hub/bw/{sample}_fw.bw", sample=SAMPLES),
    #     expand("nanoporeibat_hub/bw/{sample}_rv.bw", sample=SAMPLES)
    # params:
    #     url = "http://bioinformatik.sdu.dk/solexa/webshare/christoph/nanoporeibat_hub"
    # shell:
    #     "bin/generateUCSChub.R {input.sample_info} {params.url}"
