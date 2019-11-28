rule dexseq_prefilterIsoforms:
    input:
        annotation = "annotation/annotation.gtf",
        scaledTPM = "res/dexseq/illumina/dexseq_scaledTPM.rds",
        txdb = "annotation/annotation_txdb.sqlite",
        sampleInfo = "sample_info/sampleInfo.csv"
    params:
        threshold = 15
    output:
        "indices/dexseq/annotation_prefiltered.gtf"
    script:
        "dexseq_prefilterIsoforms.R"

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

rule dexseq_count_illumina:
    input:
        annotation = "indices/dexseq/annotation_flat.gff",
        bam = "bam/illumina/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        "dexseq/illumina/{sample}.txt"
    wildcard_constraints:
        dataSet = "illumina|teloprime"
    shell:
        "python3 bin/dexseq_count.py \
            -p yes -s reverse -f bam -r pos -a 2 \
            {input.annotation} \
            {input.bam} \
            {output}"

rule dexseq_count_nanopore:
    input:
        annotation = "indices/dexseq/annotation_flat.gff",
        bam = "bam/{dataSet}/{sample}_genome.bam"
    output:
        "dexseq/{dataSet}/{sample}.txt"
    wildcard_constraints:
        dataSet = "direct_cDNA|teloprime"
    shell:
        "python3 bin/dexseq_count.py \
            -p no -s no -f bam -r pos -a 1 \
            {input.annotation} \
            {input.bam} \
            {output}"

rule dexseq_importCounts_illumina:
    input:
        files = expand("dexseq/illumina/{sample}.txt", sample=SAMPLES),
        annotation = "indices/dexseq/annotation_flat.gff",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        dataset = "illumina"
    output:
        "res/dexseq/illumina/illumina_dxd.rds",
    script:
        "dexseq_importCounts.R"

rule dexseq_importCounts_teloprime:
    input:
        files = expand("dexseq/teloprime/{barcode}.txt", barcode=BARCODES),
        annotation = "indices/dexseq/annotation_flat.gff",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        dataset = "ont"
    output:
        "res/dexseq/teloprime/teloprime_dxd.rds",
    script:
        "dexseq_importCounts.R"

rule dexseq_diffExonUsage:
    threads:
        10
    input:
        "{file}_dxd.rds",
    output:
        dxd = "{file}_dxd_diff.rds",
        report = "{file}_dexseq.html"
    script:
        "dexseq_diffExonUsage.R"
