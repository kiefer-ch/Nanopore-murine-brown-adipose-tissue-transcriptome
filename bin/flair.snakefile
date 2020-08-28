FLAIR = config["FLAIR"]


rule flair_convert_bed12:
    input:
        bam = "bam/{dataset}/{barcode}_genome.bam",
        bai = "bam/{dataset}/{barcode}_genome.bam.bai"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    output:
        "flair/{dataset}/bed/raw/{barcode}.bed"
    shell:
        "python2 {FLAIR}/bin/bam2Bed12.py -i {input.bam} \
            > {output}"


rule flair_convert_bed12_rna:
    input:
        bam = "bam/rna/genome_{barcode}_q7_sort.bam",
        bai = "bam/rna/genome_{barcode}_q7_sort.bam.bai"
    output:
        "flair/rna/bed/raw/{barcode}.bed"
    shell:
        "python2 {FLAIR}/bin/bam2Bed12.py -i {input.bam} \
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
    elif wildcards.dataset == "rna":
        filename = "flair/rna/bed/raw/{}.bed".format(wildcards.sample)
    return filename


def get_flair_junctions(wildcards):
    if wildcards.dataset in ["cdna", "teloprime"]:
        illumina = SAMPLE_INFO.loc[wildcards.sample]["illumina"]
        filename = "flair/SJout/{}_SJ.out_filtered.tab".format(illumina)
    elif wildcards.dataset == "rna":
        if wildcards.sample == "rt":
            filename = "flair/SJout/5034_S33_SJ.out_filtered.tab"
        if wildcards.sample == "cool":
            filename = "flair/SJout/5035_S34_SJ.out_filtered.tab"
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
        dataset = "teloprime|cdna|rna"
    threads:
        4
    shell:
        "python2 {FLAIR}/flair.py correct \
            -q {input.bed} \
            -g {input.genome} \
            -j {input.junctions} \
            -c {input.chromsizes} \
            -t {threads} \
            -o {params.out_prefix}"


def get_concatNames(wildcards):
    files = list()
    if wildcards.dataset in ["teloprime", "cdna"]:
        for id in SAMPLE_INFO_ont.index:
            filename = "flair/{}/bed/corrected/{}_all_corrected.psl".format(
                wildcards.dataset, id)
            files.append(filename)
    elif wildcards.dataset == "rna":
        files = expand("flair/rna/bed/corrected/{sample}_all_corrected.psl",
            sample = ["rt", "cool"])
    return files


rule flair_concatenate:
    input:
        get_concatNames
    output:
        "flair/{dataset}/bed/corrected/concatenated_all_corrected.psl"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
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
        "TDIR=$(mktemp -d data/chip/tmp.XXXXXXXXX) && "
        "awk '{printf(\"%s\t%s\t%s\n\", $1, $2, $3)}' > $TDIR/cat.bed < {input.cage} && "
        "cat {input.combined}  >> $TDIR/cat.bed && "
        "sort -k1,1 -k2,2n $TDIR/cat.bed > $TDIR/cat.sorted.bed && "
        "bedtools merge -i $TDIR/cat.sorted.bed > {output}"


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
    elif wildcards.dataset == "rna":
        files = expand("fastq/rna/{sample}_q7.fastq.gz",
            sample = ["cool", "rt"])
    return files


rule flair_collapse:
    input:
        fastq = get_flair_fastqnames,
        genome = "annotation/genome.fa",
        annotation = "annotation/annotation.gtf",
        psl = "flair/{dataset}/bed/corrected/concatenated_all_corrected.psl",
        promoters = "data/chip/h3k4_cage_combined.bed"
    output:
        "flair/{dataset}/flair.collapse.{dataset}.isoforms.fa",
        "flair/{dataset}/flair.collapse.{dataset}.isoforms.gtf",
        "flair/{dataset}/flair.collapse.{dataset}.isoforms.psl"
    params:
        out_prefix = "flair/{dataset}/flair.collapse.{dataset}"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    threads:
        40
    shell:
        "python2 {FLAIR}/flair.py collapse \
            -g {input.genome} \
            -f {input.annotation} \
            -r {input.fastq} \
            -q {input.psl} \
            -t {threads} \
            -p {input.promoters} \
            -o {params.out_prefix} \
            -s 5 --stringent \
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
        isoforms_fasta = "flair/{dataset}/flair.collapse.{dataset}.isoforms.fa"
    output:
        "flair/{dataset}/flair_{dataset}_counts_matrix.tsv"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    threads:
        40
    shell:
        "python2 {FLAIR}/flair.py quantify \
            -r {input.reads_manifest} \
            -i {input.isoforms_fasta} \
            -o {output} \
            -t {threads} \
            --temp_dir ./ \
            --trust_ends"


rule flair_example_browser_track:
    input:
        txdb = "annotation/annotation_txdb.sqlite",
        txdb_flair = "flair/cdna/flair.collapse.cdna.isoforms_txdb.sqlite",
        illumina_warm = "bam/illumina/5034_S33_Aligned.sortedByCoord.out.bam",
        illumina_cold = "bam/illumina/5035_S34_Aligned.sortedByCoord.out.bam",
        cdna_warm = "bam/cdna/barcode07_genome.bam",
        cdna_cold = "bam/cdna/barcode08_genome.bam",
        genome = "annotation/genome.fa",
    output:
        "res/comparisons/flair_browser_examples.html"
    script:
        "flair_example_browser_track.Rmd"
