# salmon align_reannotation
BEDTOOLS = config["BEDTOOLS"]
MASHMAP = config["MASHMAP"]


rule salmon_prepareDecoys_reannotatedStringtie:
    input:
        genome = "data/annotation/genome.fa",
        annotation = "data/reannotation/stringtie/{dataset}_stringtie_nophuchatooStrand.gtf",
        transcripts = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.fa"
    output:
        "indices/salmon_{dataset}_stringtie/decoy/gentrome.fa",
        "indices/salmon_{dataset}_stringtie/decoy/decoys.txt"
    params:
        outputDir = "indices/salmon_{dataset}_stringtie/decoy"
    threads: 20
    wildcard_constraints:
        dataset = "illumina|teloprime"
    shell:
        "bin/generateDecoyTranscriptome.sh \
            -j {threads} \
            -g {input.genome} \
            -t {input.transcripts} \
            -a {input.annotation} \
            -b {BEDTOOLS} \
            -m {MASHMAP} \
            -o {params.outputDir}"


rule salmon_prepareDecoys_reannotatedFlair:
    input:
        genome = "data/annotation/genome.fa",
        annotation = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.gtf",
        transcripts = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.fa"
    output:
        temp("indices/salmon_{dataset}_flair/decoy/gentrome.fa"),
        temp("indices/salmon_{dataset}_flair/decoy/decoys.txt")
    params:
        outputDir = "indices/salmon_{dataset}_flair/decoy"
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


rule salmon_index_reannotated:
    input:
        gentrome = "indices/salmon_{dataset}_{method}/decoy/gentrome.fa",
        decoy = "indices/salmon_{dataset}_{method}/decoy/decoys.txt"
    output:
        "indices/salmon_{dataset}_{method}/duplicate_clusters.tsv",
        "indices/salmon_{dataset}_{method}/hash.bin",
        "indices/salmon_{dataset}_{method}/rsd.bin",
        "indices/salmon_{dataset}_{method}/sa.bin",
        "indices/salmon_{dataset}_{method}/txpInfo.bin"
    params:
        outputDir = "indices/salmon_{dataset}_{method}"
    wildcard_constraints:
        method = "flair|stringtie",
        dataset = "illumina|cdna|teloprime"
    threads: 20
    shell:
        "salmon index \
            --gencode \
            -t {input.gentrome} \
            -i {params.outputDir} \
            -d {input.decoy} \
            -p {threads}"


rule salmon_align_reannotated:
    threads: 5
    input:
        "indices/salmon_{dataset}_{method}/duplicate_clusters.tsv",
        "indices/salmon_{dataset}_{method}/hash.bin",
        "indices/salmon_{dataset}_{method}/rsd.bin",
        "indices/salmon_{dataset}_{method}/sa.bin",
        "indices/salmon_{dataset}_{method}/txpInfo.bin",
        fastq_fw = "data/fastq/illumina/trimmed/{sample}_R1_001_trimmed.fastq.gz",
        fastq_rv = "data/fastq/illumina/trimmed/{sample}_R2_001_trimmed.fastq.gz"
    output:
        "data/quantification/illumina_{dataset}_{method}/{sample}/quant.sf"
    params:
        inputDir = "indices/salmon_{dataset}_{method}",
        outputDir = "data/quantification/illumina_{dataset}_{method}/{sample}"
    wildcard_constraints:
        method = "flair|stringtie"
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


# minimap
def get_minimapIndexInput(wildcards):
    if wildcards.method == "flair":
        file = "data/reannotation/flair/annotation/{dataset}_flair.isoforms.fa"
    elif wildcards.method == "stringtie":
        file = "data/reannotation/stringtie/{dataset}_stringtie_noUnknownStrand.fa"
    return file


rule minimap_index_reannotation:
    input:
        get_minimapIndexInput
    output:
        "indices/minimap2_{dataset}_{method}/minimap2.mmi"
    threads:
        3
    conda:
        "../envs/nanopore.yaml"
    shell:
        "minimap2 \
        -t {threads} \
        -d {output} \
        {input}"


def get_minimapInput(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "data/fastq/{}/{}/{}_q7.fastq.gz".format(
            wildcards.dataset, wildcards.flowcell, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "data/fastq/rna/{}_q7.fastq.gz".format(wildcards.barcode)
    return file_name


rule minimap_mapReannotation:
    input:
        index = "indices/minimap2_{dataset}_{method}/minimap2.mmi",
        fastq = "data/fastq/cdna/{flowcell}/{barcode}_q7.fastq.gz"
    output:
        sam = temp("data/bam/cdna_{dataset}_{method}/{flowcell}_{barcode}_sort.bam")
    threads:
        20
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        type = "flair|stringtie"
    conda:
        "../envs/nanopore.yaml"
    shell:
        "minimap2 \
            -2 \
            -ax map-ont \
            -t {threads} \
            --secondary=no \
            -uf \
            {input.index} \
            {input.fastq} | \
        samtools sort -l 5 -o {output} -O bam -@ 6"


rule samtools_merge2:
    input:
        "data/bam/cdna_{dataset}_{method}/flowcell1_{barcode}_sort.bam",
        "data/bam/cdna_{dataset}_{method}/flowcell2_{barcode}_sort.bam"
    output:
        temp("data/reannotation/cdna_{dataset}_{method}/{barcode}_sort.bam")
    threads:
        6
    conda:
        "../envs/nanopore.yaml"
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        method = "stringtie|flair",
        barcode = "barcode.."
    shell:
        "samtools merge -@ {threads} {output} {input}"


rule quantify_minimap_reannotated:
    input:
        "data/reannotation/cdna_{dataset}_{method}/{barcode}_sort.bam"
    output:
        "data/quantification/cdna_{dataset}_{method}/{barcode}_quant.tsv"
    wildcard_constraints:
        dataset = "teloprime|cdna|illumina",
        type = "flair|stringtie"
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "quantify_minimap.R"
