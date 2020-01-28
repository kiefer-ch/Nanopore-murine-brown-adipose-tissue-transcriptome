rule tximport_dexseq_illumina:
    input:
        salmon_out = expand("salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/dexseq/illumina/dexseq_scaledTPM.rds"
    script:
        "dexseq_txImport.R"


rule dexseq_prefilterIsoforms_illumina:
    input:
        annotation = "annotation/annotation.gtf",
        scaledTPM = "res/dexseq/illumina/dexseq_scaledTPM.rds",
        txdb = "annotation/annotation_txdb.sqlite",
        sampleInfo = "sample_info/sampleInfo.csv"
    params:
        threshold = 15
    output:
        "indices/dexseq/annotation_illumina_prefiltered.gtf"
    script:
        "dexseq_prefilterIsoforms.R"


rule dexseq_prefilterIsoforms_ont:
    input:
        counts = "res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_ntd.csv.gz",
        annotation = "annotation/annotation.gtf",
        txdb = "annotation/annotation_txdb.sqlite",
        sampleInfo = "sample_info/sampleInfo.csv"
    params:
        threshold = 15
    output:
        "indices/dexseq/annotation_{dataset}_prefiltered.gtf"
    wildcard_constraints:
        dataset = "teloprime|cdna"
    script:
        "dexseq_prefilterIsoforms_ont.R"


rule dexseq_prepareAnnotation:
    input:
        "indices/dexseq/annotation_{dataset}_prefiltered.gtf"
    output:
        gff = "indices/dexseq/annotation_{dataset}_flat.gff",
        gtf = "indices/dexseq/annotation_{dataset}_flat.gtf"
    wildcard_constraints:
        dataset = "illumina|teloprime|cdna"
    shell:
        "python3 bin/dexseq_prepare_annotation2.py \
            --aggregate no \
            -f {output.gtf}\
            {input} \
            {output.gff}"


def get_bam_ont(wildcards):
    files = list()
    if wildcards.dataset == "teloprime":
        for barcode in SAMPLE_INFO_ont["ont"]:
            filename = "bam/{}/{}_genome.bam".format(
                wildcards.dataset, barcode)
            files.append(filename)
    elif wildcards.dataset == "cdna":
        for barcode in SAMPLE_INFO_ont["cdna"]:
            filename = "bam/{}/{}_genome.bam".format(
                wildcards.dataset, barcode)
            files.append(filename)
    return files


rule featureCounts_count_ont:
    input:
        files = get_bam_ont,
        annotation = "indices/dexseq/annotation_{dataset}_flat.gtf",
    threads:
        20
    output:
        "res/dexseq/{dataset}/{dataset}_featureCounts.out"
    wildcard_constraints:
        dataset = "teloprime|cdna"
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
                       sample=SAMPLES_ont),
        annotation = "indices/dexseq/annotation_illumina_flat.gtf",
    threads:
        20
    output:
        "res/dexseq/illumina/illumina_featureCounts.out"
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
        counts = "res/dexseq/{dataset}/{dataset}_featureCounts.out",
        annotation = "indices/dexseq/annotation_{dataset}_flat.gtf",
        sample_info = "sample_info/sampleInfo.csv"
    wildcard_constraints:
        dataset = "illumina|teloprime|cdna"
    output:
        "res/dexseq/{dataset}/{dataset}_dxd.rds",
    script:
        "dexseq_importCounts_from_featureCounts.R"


rule dexseq_diffExonUsage:
    threads:
        10
    input:
        "res/dexseq/{dataset}/{dataset}_dxd.rds",
    output:
        dxd = "res/dexseq/{dataset}/{dataset}_dxd_diff.rds",
        report = "res/dexseq/{dataset}/{dataset}_dexseq.html"
    wildcard_constraints:
        dataset = "illumina|teloprime|cdna"
    script:
        "dexseq_diffExonUsage.R"


rule dexseq_heatmap:
    input:
        dxd = "res/dexseq/{dataset}/{dataset}_dxd_diff.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    params:
        pvalue_cutoff = .05
    output:
        "res/dexseq/{dataset}/{dataset}_heatmap.html"
    wildcard_constraints:
        dataset = "illumina|teloprime|cdna"
    script:
        "dexseq_heatmap.Rmd"

rule dexseq_resultsTable:
    input:
        dxd = "res/dexseq/{dataset}/{dataset}_dxd_diff.rds",
        biomaRt_gene = "annotation/biomaRt_gene.rds"
    output:
        "res/dexseq/{dataset}/{dataset}_dexseq_results.csv.gz"
    wildcard_constraints:
        dataset = "illumina|teloprime|cdna"
    script:
        "dexseq_resultsTable.R"
