rule deseq_exportCmONT:
    input:
        dds = "data/deseq/{dataset}/{dataset}_dds_{type}.rds"
    params:
        type = "{type}"
    wildcard_constraints:
        level = "gene|transcript",
        dataset = "rna|cdna|teloprime"
    output:
        cts = "data/deseq/{dataset}/{dataset}_{type}_cm_cts.csv.gz",
        ntd = "data/deseq/{dataset}/{dataset}_{type}_cm_ntd.csv.gz"
    script:
        "deseq_exportCm.R"


rule deseq_exportCmIllumina:
    input:
        dds = "data/deseq/{dataset}/{dataset}_dds_{type}.rds"
    params:
        type = "{type}"
    wildcard_constraints:
        level = "gene|transcript",
        dataset = "illumina"
    output:
        cts = "data/deseq/{dataset}/{dataset}_{type}_cm_cts.csv.gz",
        ntd = "data/deseq/{dataset}/{dataset}_{type}_cm_ntd.csv.gz",
        tpm = "data/deseq/{dataset}/{dataset}_{type}_cm_tpm.csv.gz"
    script:
        "deseq_exportCm.R"



# DGE
def get_biomart(wildcards):
    if wildcards.level == "genelevel":
        return "annotation/biomaRt_gene.rds"
    elif wildcards.level == "txlevel":
        return "annotation/biomaRt_tx.rds"


def get_threshold(wildcards):
    if wildcards.shrink == "apeglm":
        return 1
    elif wildcards.shrink == "noShrink":
        return 0

rule deseq_dge:
    input:
        dds = "{file_path}/{dataset}_{level}_dds.rds",
        biomart = get_biomart
    output:
        "{file_path}/{dataset}_{level}_{shrink}_results.csv.gz",
    params:
        level = "{level}",
        shrink = "{shrink}",
        lfcThreshold = get_threshold,
        alpha = 0.05
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    threads: 4
    script:
        "deseq_dge.R"


# pathway analysis
rule topgo_analysis:
    input:
        "{file}_{level}_{shrink}_results.csv.gz"
    output:
        "{file}_{level}_{shrink}_topgo.rds"
    params:
        level = "{level}",
        shrink = "{shrink}"
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    script:
        "deseq_topgo.R"


rule reactome_analysis:
    input:
        "{file}_{level}_{shrink}_results.csv.gz"
    output:
        "{file}_{level}_{shrink}_reactome.rds",
        "{file}_{level}_{shrink}_reactomeGL.rds"
    params:
        level = "{level}",
        shrink = "{shrink}"
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    script:
        "deseq_reactome.R"


# report
rule deseq_report:
    input:
        dds = "{file}_{level}_dds.rds",
        rld = "{file}_{level}_cm_rld.csv.gz",
        res = "{file}_{level}_{shrink}_results.csv.gz",
        topgo = "{file}_{level}_{shrink}_topgo.rds",
        reactome = "{file}_{level}_{shrink}_reactome.rds",
        reactome_genelist = "{file}_{level}_{shrink}_reactomeGL.rds",
        grouped_biotypes = "data/biotype_groups.csv"
    params:
        level = "{level}",
        shrink = "{shrink}",
        lfcThreshold = get_threshold,
        alpha = .05
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    output:
        "{file}_{level}_{shrink}_report.html"
    script:
        "deseq_report.Rmd"
