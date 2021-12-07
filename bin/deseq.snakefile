# count matrix export
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


rule deseq_dge:
    input:
        dds = "data/deseq/{dataset}/{dataset}_dds_{type}.rds",
        biomart = get_biomart
    output:
        "data/deseq/{dataset}/{dataset}_deseqResults_{type}.tsv.gz",
    params:
        level = "{type}",
        lfcThreshold = 0.5,
        alpha = 0.05
    wildcard_constraints:
        level = "genelevel|txlevel",
        shrink = "apeglm|noShrink"
    threads: 4
    script:
        "deseq_dge.R"
