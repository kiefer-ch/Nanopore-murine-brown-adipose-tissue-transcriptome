rule flair_convert_bed12:
    input:
        "bam/{dataset}/{barcode}_genome.bam"
    wildcard_constraints:
        dataset = "teloprime"
    output:
        "flair/{dataset}/bed/raw/{barcode}.bed"
    shell:
        "python2 ~/src/flair/bin/bam2Bed12.py -i {input} \
            > {output}"


rule filter_SJout:
    input:
        "bam/illumina/{sample}_SJ.out.tab"
    output:
        "flair/SJout/{sample}_SJ.out_filtered.tab"
    params:
        threshold = 10
    script:
        "filter_SJout.R"


rule flair_correct:
    input:
        bed = "flair/{dataset}/bed/raw/{barcode}.bed",
        genome = "annotation/genome.fa",
        chromsizes = "annotation/genome.fa.fai",
        junctions = "flair/SJout/{sample}_SJ.out_filtered.tab"
    output:
        "flair/{dataset}/bed/corrected/{sample}_{barcode}_all_corrected.bed",
        "flair/{dataset}/bed/corrected/{sample}_{barcode}_all_corrected.psl",
        "flair/{dataset}/bed/corrected/{sample}_{barcode}_all_inconsistent.bed"
    params:
        out_prefix = "flair/{dataset}/bed/corrected/{sample}_{barcode}"
    wildcard_constraints:
        dataset = "teloprime",
        barcode = "barcode0[1-6]",
        sample = "50[34][0-9]_S[3-4][0-9]"
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


def get_flair_filenames_teloprime():
    files = list()
    for i in range(0, len(SAMPLES_ont)):
        filename = "flair/teloprime/bed/corrected/{}_{}_all_corrected.psl".format(SAMPLES_ont[i], BARCODES[i])
        files.append(filename)
    return files


rule flair_concatenate_teloprime:
    input:
        get_flair_filenames_teloprime()
    output:
        "flair/teloprime/bed/corrected/concatenated_all_corrected.psl"
    shell:
        "cat {input} > {output}"


rule flair_collapse:
    input:
        genome = "annotation/genome.fa",
        annotation = "annotation/annotation.gtf",
        psl = "flair/{dataset}/bed/corrected/concatenated_all_corrected.psl",
        fastq = expand("fastq/teloprime/{flowcell}/{barcode}_q7.fastq.gz",
                       flowcell=["X1_flowcell", "X3_flowcell"], barcode=BARCODES)
    output:
        "flair/{dataset}/flair.collapse.isoforms.fa",
        "flair/{dataset}/flair.collapse.isoforms.gtf",
        "flair/{dataset}/flair.collapse.isoforms.psl"
    params:
        out_prefix = "flair/{dataset}/flair.collapse"
    wildcard_constraints:
        dataset = "teloprime"
    threads:
        40
    shell:
        "python2 ~/src/flair/flair.py collapse \
            -g {input.genome} \
            -f {input.annotation} \
            -r {input.fastq} \
            -q {input.psl} \
            -t {threads} \
            -o {params.out_prefix} \
            -s 10 \
            --temp_dir ./"
