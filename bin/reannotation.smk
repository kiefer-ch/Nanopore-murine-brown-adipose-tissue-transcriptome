FLAIR = config["FLAIR"]


rule flair_convert_bed12:
    input:
        "data/bam/{dataset}/merged/{dataset}_{barcode}_genome_primaryOnly.bam.bai",
        bam = "data/bam/{dataset}/merged/{dataset}_{barcode}_genome_primaryOnly.bam"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    output:
        temp("data/reannotation/flair/bed/raw/{dataset}/{dataset}_{barcode}.bed")
    shell:
        "python2 {FLAIR}/bin/bam2Bed12.py -i {input.bam} \
            > {output}"


rule filter_SJout:
    input:
        "data/bam/illumina/{sample}_SJ.out.tab"
    output:
        temp("data/reannotation/flair/SJout/{sample}_SJ.out_filtered.tab")
    params:
        threshold = config["SJ_cutoff"]
    script:
        "reannotation_filterSJout.R"


def get_flair_junctions(wildcards):
    if wildcards.dataset == "cdna":
        illumina = SAMPLE_INFO_ont.set_index("cdna").loc[wildcards.barcode, "illumina"]
        filename = "data/reannotation/flair/SJout/{}_SJ.out_filtered.tab".format(illumina)
    elif wildcards.dataset == "teloprime":
        illumina = SAMPLE_INFO_ont.set_index("ont").loc[wildcards.barcode, "illumina"]
        filename = "data/reannotation/flair/SJout/{}_SJ.out_filtered.tab".format(illumina)
    elif wildcards.dataset == "rna":
        if wildcards.barcode == "rt":
            filename = "data/reannotation/flair/SJout/5034_S33_SJ.out_filtered.tab"
        if wildcards.barcode == "cool":
            filename = "data/reannotation/flair/SJout/5035_S34_SJ.out_filtered.tab"
    return filename


rule flair_correct:
    input:
        bed = "data/reannotation/flair/bed/raw/{dataset}/{dataset}_{barcode}.bed",
        junctions = get_flair_junctions,
        genome = "data/annotation/genome.fa",
        chromsizes = "data/annotation/genome.fa.fai"
    output:
        temp(multiext("data/reannotation/flair/bed/corrected/{dataset}/{dataset}_{barcode}_all_corrected",
            ".bed", ".psl")),
        temp("data/reannotation/flair/bed/corrected/{dataset}/{dataset}_{barcode}_all_inconsistent.bed")
    params:
        out_prefix = "data/reannotation/flair/bed/corrected/{dataset}/{dataset}_{barcode}"
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
    if wildcards.dataset == "teloprime":
        for id in SAMPLE_INFO_ont["ont"]:
            filename = "data/reannotation/flair/bed/corrected/{}/{}_{}_all_corrected.psl".format(
                wildcards.dataset, wildcards.dataset, id)
            files.append(filename)
    elif wildcards.dataset == "cdna":
        for id in SAMPLE_INFO_ont["cdna"]:
            filename = "data/reannotation/flair/bed/corrected/{}/{}_{}_all_corrected.psl".format(
                wildcards.dataset, wildcards.dataset, id)
            files.append(filename)
    elif wildcards.dataset == "rna":
        files = expand("data/reannotation/flair/bed/corrected/rna/rna_{sample}_all_corrected.psl",
            sample = ["rt", "cool"])
    return files


rule flair_concatenate:
    input:
        get_concatNames
    output:
        "data/reannotation/flair/bed/corrected/{dataset}/{dataset}_concatenated_all_corrected.psl"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    shell:
        "cat {input} > {output}"


rule bedtools_combineH3k4me3:
    input:
        expand("data/chip/k4me3/{condition}_peaks.narrowPeak",
            condition=["cold_hfd", "cold_ncd", "warm_hfd", "warm_ncd"])
    output:
        "data/chip/k4me3/combined_peaks.narrowPeak"
    shell:
        """
        TDIR=$(mktemp -d data/chip/k4me3/tmp.XXXXXXXXX) &&
        cat {input} > $TDIR/cat.bed &&
        sort -k1,1 -k2,2n $TDIR/cat.bed > $TDIR/cat.sorted.bed &&
        bedtools merge -i $TDIR/cat.sorted.bed > {output}
        rm -r $TDIR
        """


rule bedtools_combineH3k4me3Cage:
    input:
        combined = "data/chip/k4me3/combined_peaks.narrowPeak",
        # from FANTOM
        cage = "data/chip/cage/mm10.cage_peak_phase1and2combined_coord.bed"
    output:
        "data/chip/h3k4_cage_combined.bed"
    shell:
        """
        TDIR=$(mktemp -d data/chip/tmp.XXXXXXXXX) &&
        awk '{{printf(\"%s\\t%s\\t%s\\n\", $1, $2, $3)}}' > $TDIR/cat.bed < {input.cage} &&
        cat {input.combined} >> $TDIR/cat.bed &&
        sort -k1,1 -k2,2n $TDIR/cat.bed > $TDIR/cat.sorted.bed &&
        bedtools merge -i $TDIR/cat.sorted.bed > {output}
        rm -r $TDIR
        """


def get_flair_fastqnames(wildcards):
    files = list()
    if wildcards.dataset == "teloprime":
        for barcode in SAMPLE_INFO_ont["ont"]:
            for flowcell in ["flowcell1", "flowcell2"]:
                filename = "data/fastq/{}/{}/{}_q7.fastq.gz".format(
                    wildcards.dataset, flowcell, barcode)
                files.append(filename)
    elif wildcards.dataset == "cdna":
        for barcode in SAMPLE_INFO_ont["cdna"]:
            for flowcell in ["flowcell1", "flowcell2"]:
                filename = "data/fastq/{}/{}/{}_q7.fastq.gz".format(
                    wildcards.dataset, flowcell, barcode)
                files.append(filename)
    elif wildcards.dataset == "rna":
        files = expand("data/fastq/rna/{sample}_q7.fastq.gz",
            sample = ["cool", "rt"])
    return files


def ends_switch(wildcards):
    par = ""
    if wildcards.dataset == "teloprime":
        par = "--trust_ends"
    elif wildcards.dataset == "cdna":
        par = ""
    return par


rule flair_collapse:
    input:
        fastq = get_flair_fastqnames,
        genome = "data/annotation/genome.fa",
        annotation = "data/annotation/annotation.gtf",
        psl = "data/reannotation/flair/bed/corrected/{dataset}/{dataset}_concatenated_all_corrected.psl",
        promoters = "data/chip/h3k4_cage_combined.bed"
    output:
        multiext("data/reannotation/flair/annotation/{dataset}_flair.isoforms",
            ".fa", ".gtf", ".psl")
    params:
        out_prefix = "data/reannotation/flair/annotation/{dataset}_flair",
        ends = ends_switch
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
            {params.ends} \
            --temp_dir ./"















rule flair_quantify:
    input:
        get_flair_fastqnames,
        reads_manifest = "sample_info/flair_{dataset}_readsManifest.tsv",
        isoforms_fasta = "flair/{dataset}/flair.collapse.{dataset}.isoforms.fa"
    output:
        "flair/{dataset}/flair_{dataset}_counts_matrix.tsv"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    params:
        ends = ends_switch
    threads:
        40
    shell:
        "python2 {FLAIR}/flair.py quantify \
            -r {input.reads_manifest} \
            -i {input.isoforms_fasta} \
            -o {output} \
            -t {threads} \
            {params.ends} \
            --temp_dir ./"
