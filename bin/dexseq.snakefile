rule dexseq_prefilterIsoforms:
    input:
        annotation = "annotation/annotation.gtf",
        scaledTPM = "data/scaledTPM_all.rds",
        txdb = "annotation/annotation.txdb.rds",
        sampleInfo = "sample_info/sampleInfo.csv"
    params:
        threshold = 15
    output:
        "indices/dexseq/annotation_prefiltered.gtf"
    script:
        "dexseq_prefilterIsoforms.R"

# https://github.com/vivekbhr/Subread_to_DEXSeq
rule dexseq_prepareAnnotation:
    input:
        "indices/dexseq/annotation_prefiltered.gtf"
    output:
        "indices/dexseq/annotation_flat.gff"
    shell:
        "python3 bin/dexseq_prepare_annotation.py \
            --aggregate no \
            {input} \
            {output}"

rule dexseq_count:
    input:
        annotation = "indices/dexseq/annotation_flat.gff",
        bam = "BAM/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "dexseq/{sample}.txt"
    shell:
        "python3 bin/dexseq_count.py \
            -p yes -s yes -f bam -r pos \
            {input.annotation} \
            {input.bam} \
            {output}"

rule dexseq_importCounts:
    input:
        expand("dexseq/{sample}.txt", sample=SAMPLES),
        annotation = "indices/dexseq/annotation_flat.gff",
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "data/dxd.rds",
    script:
        "dexseq_importCounts.R"

rule dexseq_analyse:
    threads:
        4
    input:
        "data/dxd.rds"
    output:
        "res/dexseq/desxseq.html"
    script:
        "dexseq_analysis.Rmd"
