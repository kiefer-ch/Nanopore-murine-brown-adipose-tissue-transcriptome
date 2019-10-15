# https://github.com/vivekbhr/Subread_to_DEXSeq
rule dexseq_prepare_annotation:
    input:
        "annotation/annotation.gtf"
    output:
        "indices/dexseq/annotation_flat.gff"
    shell:
        "python3 bin/dexseq_prepare_annotation.py \
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
#        expand("dexseq/{sample}.txt", sample=SAMPLES),
        annotation = "indices/dexseq/annotation_flat.gff",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        filter_by = "ont"
    output:
        "data/dds_gencode.vM22_gene_ont.rds"
    script:
        "dexseq_importCounts.R"
