rule dexseq_prefilterIsoforms:
    input:
        annotation = "annotation/annotation.gtf",
        scaledTPM = "data/scaledTPM_all.rds"
    threads:
        2
    output:
        "desxseq/annotation_prefiltered.gtf"
    script:
        "dexseq_prefilterIsoforms.R"

# https://github.com/vivekbhr/Subread_to_DEXSeq
rule dexseq_prepare_annotation:
    input:
        "annotation/annotation.gtf"
    output:
        "indices/dexseq/annotation_flat.gff"
    shell:
        "python3 bin/dexseq_prepare_annotation.py \
            --aggregate no \
            {input} \
            {output}"

rule dexseq_count_illumina:
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
