rule flair_convert_bed12:
    input:
        bam = "bam/{dataset}/{barcode}_genome.bam",
        bai = "bam/{dataset}/{barcode}_genome.bam.bai"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    output:
        "flair/{dataset}/bed/raw/{barcode}.bed"
    shell:
        "python2 ~/src/flair/bin/bam2Bed12.py -i {input.bam} \
            > {output}"


rule filter_SJout:
    input:
        "bam/illumina/{sample}_SJ.out.tab"
    output:
        "flair/SJout/{sample}_SJ.out_filtered.tab"
    params:
        threshold = 10
    script:
        "flair_filterSJout.R"


def get_flair_bednames(wildcards):
    if wildcards.dataset == "cdna":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["cdna"]
        filename = "flair/cdna/bed/raw/{}.bed".format(barcode)
    elif wildcards.dataset == "teloprime":
        barcode = SAMPLE_INFO.loc[wildcards.sample]["ont"]
        filename = "flair/teloprime/bed/raw/{}.bed".format(barcode)
    return filename


def get_flair_junctions(wildcards):
    illumina = SAMPLE_INFO.loc[wildcards.sample]["illumina"]
    filename = "flair/SJout/{}_SJ.out_filtered.tab".format(illumina)
    return filename


rule flair_correct:
    input:
        bed = get_flair_bednames,
        junctions = get_flair_junctions,
        genome = "annotation/genome.fa",
        chromsizes = "annotation/genome.fa.fai"
    output:
        "flair/{dataset}/bed/corrected/{sample}_all_corrected.bed",
        "flair/{dataset}/bed/corrected/{sample}_all_corrected.psl",
        "flair/{dataset}/bed/corrected/{sample}_all_inconsistent.bed"
    params:
        out_prefix = "flair/{dataset}/bed/corrected/{sample}"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    threads:
        4
    shell:
        "python2 ~/src/flair/flair.py correct \
            -q {input.bed} \
            -g {input.genome} \
            -j {input.junctions} \
            -c {input.chromsizes} \
            -t {threads} \
            -o {params.out_prefix}"


def get_concatNames(wildcards):
    files = list()
    for id in SAMPLE_INFO_ont.index:
        filename = "flair/{}/bed/corrected/{}_all_corrected.psl".format(
            wildcards.dataset, id)
        files.append(filename)
    return files


rule flair_concatenate:
    input:
        get_concatNames
    output:
        "flair/{dataset}/bed/corrected/concatenated_all_corrected.psl"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    shell:
        "cat {input} > {output}"


rule bedtools_combineWarmColdH3k4:
    input:
        "data/chip/k4me3/NC_broad.bed"
        "data/chip/k4me3/NW_broad.bed"
    output:
        "data/chip/k4me3/combined_broad.bed"
    shell:
        "TDIR=$(mktemp -d data/chip/k4me3/tmp.XXXXXXXXX) \
        cat {input} > $TDIR/cat.bed \
        sort -k1,1 -k2,2n $TDIR/cat.bed > $TDIR/cat.sorted.bed \
        bedtools merge -i $TDIR/cat.sorted.bed > {output}"


rule bedtools_combine_h3k4me3_cage:
    input:
        combined = "data/chip/k4me3/combined_broad.bed",
        # from FANTOM
        cage = "data/chip/cage/mm10.cage_peak_phase1and2combined_coord.bed"
    output:
        "data/chip/h3k4_cage_combined.bed"
    shell:
        "TDIR=$(mktemp -d data/chip/tmp.XXXXXXXXX) && \
        awk '{printf ("%s\t%s\t%s\n", $1, $2, $3)}' > $TDIR/cat.bed < {input.cage}
        cat {input.combined}  >> $TDIR/cat.bed && \
        sort -k1,1 -k2,2n $TDIR/cat.bed > $TDIR/cat.sorted.bed && \
        bedtools merge -i $TDIR/cat.sorted.bed > {output}"

def get_flair_fastqnames(wildcards):
    files = list()
    if wildcards.dataset == "teloprime":
        for barcode in SAMPLE_INFO_ont["ont"]:
            filename = "fastq/{}/merged/{}_merged.fastq.gz".format(
                wildcards.dataset, barcode)
            files.append(filename)
    elif wildcards.dataset == "cdna":
        for barcode in SAMPLE_INFO_ont["cdna"]:
            filename = "fastq/{}/merged/{}_merged.fastq.gz".format(
                wildcards.dataset, barcode)
            files.append(filename)
    return files


rule flair_collapse:
    input:
        fastq = get_flair_fastqnames,
        genome = "annotation/genome.fa",
        annotation = "annotation/annotation.gtf",
        psl = "flair/{dataset}/bed/corrected/concatenated_all_corrected.psl",
        promoters = "data/chip/k4me3/combined_broad.bed"
    output:
        "flair/{dataset}/flair.collapse.isoforms.fa",
        "flair/{dataset}/flair.collapse.isoforms.gtf",
        "flair/{dataset}/flair.collapse.isoforms.psl"
    params:
        out_prefix = "flair/{dataset}/flair.collapse"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    threads:
        40
    shell:
        "python2 ~/src/flair/flair.py collapse \
            -g {input.genome} \
            -f {input.annotation} \
            -r {input.fastq} \
            -q {input.psl} \
            -t {threads} \
            -p {input.promoters} \
            -o {params.out_prefix} \
            -s 10 --stringent \
            --temp_dir ./"


rule merge_fastq_teloprime:
    input:
        X1 = "fastq/teloprime/X1_flowcell/{barcode}_q7.fastq.gz",
        X3 = "fastq/teloprime/X3_flowcell/{barcode}_q7.fastq.gz"
    output:
        "fastq/teloprime/merged/{barcode}_merged.fastq.gz"
    shell:
        "cat {input.X1} {input.X3} > {output}"


rule merge_fastq_cdna:
    input:
        X1 = "fastq/cdna/pool1/{barcode}_q7.fastq.gz",
        X3 = "fastq/cdna/pool2/{barcode}_q7.fastq.gz"
    output:
        "fastq/cdna/merged/{barcode}_merged.fastq.gz"
    shell:
        "cat {input.X1} {input.X3} > {output}"


rule flair_quantify:
    input:
        get_flair_fastqnames,
        reads_manifest = "sample_info/flair_{dataset}_readsManifest.tsv",
        isoforms_fasta = "flair/{dataset}/flair.collapse.isoforms.fa"
    output:
        "flair/{dataset}/flair_{dataset}_counts_matrix.tsv"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    threads:
        40
    shell:
        "python2 ~/src/flair/flair.py quantify \
            -r {input.reads_manifest} \
            -i {input.isoforms_fasta} \
            -o {output} \
            -t {threads} \
            --temp_dir ./ \
            --trust_ends"
