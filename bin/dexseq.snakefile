rule dexseq_prefilterIsoforms_illumina:
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


rule dexseq_prefilterIsoforms_teloprime:
    input:
        annotation = "annotation/annotation.gtf",
        txdb = "annotation/annotation_txdb.sqlite",
        counts = "res/deseq/teloprime/txlevel/teloprime_txlevel_cm_ntd.csv.gz",
                sampleInfo = "sample_info/sampleInfo.csv"
    params:
        threshold = 15
    output:
        "indices/dexseq/annotation_teloprime_prefiltered.gtf"
    script:
        "dexseq_prefilterIsoforms_ont.R"


rule dexseq_prepareAnnotation:
    input:
        "indices/dexseq/annotation_{dataset}_prefiltered.gtf"
    output:
        gff = "indices/dexseq/annotation_{dataset}_flat.gff",
        gtf = "indices/dexseq/annotation_{dataset}_flat.gtf"
    wildcard_constraints:
        dataset = "illumina|teloprime"
    shell:
        "python3 bin/dexseq_prepare_annotation2.py \
            --aggregate no \
            -f {output.gtf}\
            {input} \
            {output.gff}"


rule featureCounts_count_teloprime:
    input:
        files = expand("bam/teloprime/{barcode}_genome.bam", barcode=BARCODES),
        annotation = "indices/dexseq/annotation_teloprime_flat.gtf",
    threads:
        20
    output:
        "res/dexseq/teloprime/teloprime_featureCounts.out"
    shell:
        "featureCounts --donotsort \
            -L \
            -f \
            -O \
            -s 0 \
            -T {threads} \
            -F GTF \
            -a {input.annotation} \
            -o {output} \
            {input.files}"

rule featureCounts_count_teloprime_flair:
    input:
        files = expand("bam/teloprime/{barcode}_genome.bam", barcode=BARCODES),
        annotation = "flair/teloprime/flair.collapse.isoforms.gtf",
    threads:
        20
    output:
        "res/dexseq/teloprime_flair/teloprime_flair_featureCounts.out"
    shell:
        "featureCounts --donotsort \
            -L \
            -f \
            -O \
            -s 0 \
            -T {threads} \
            -F GTF \
            -a {input.annotation} \
            -o {output} \
            {input.files}"


rule featureCounts_count_illumina:
    input:
        files = expand("bam/illumina/{sample}_Aligned.sortedByCoord.out.bam",
                       sample=SAMPLES),
        annotation = "indices/dexseq/annotation_illumina_flat.gtf",
    threads:
        20
    output:
        "res/dexseq/illumina/illlumina_featureCounts.out"
    shell:
        "featureCounts --donotsort \
            -s 2 -p \
            -f \
            -O \
            -T {threads} \
            -F GTF \
            -a {input.annotation} \
            -o {output} \
            {input.files}"


rule dexseq_importCounts_from_featureCounts:
    input:
        counts = "res/dexseq/{dataset}/{dataset}_flair_featureCounts.out",
        annotation = "indices/dexseq/annotation_{dataset}_flat.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    wildcard_constraints:
        dataset = "illumina|teloprime"
    output:
        "res/dexseq/{dataset}/{dataset}_dxd.rds",
    script:
        "dexseq_importCounts_from_featureCounts.R"

rule dexseq_importCounts_from_featureCounts_flair:
    input:
        counts = "res/dexseq/{dataset}_flair/{dataset}_featureCounts.out",
        annotation = "flair/{dataset}/flair.collapse.isoforms.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    wildcard_constraints:
        dataset = "illumina|teloprime"
    output:
        "res/dexseq/{dataset}_flair/{dataset}_flair_dxd.rds",
    script:
        "dexseq_importCounts_from_featureCounts.R"

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

rule dexseq_heatmap:
    input:
        dxd = "{file}_dxd_diff.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    params:
        pvalue_cutoff = .05
    output:
        "{file}_heatmap.html"
    script:
        "dexseq_heatmap.Rmd"
