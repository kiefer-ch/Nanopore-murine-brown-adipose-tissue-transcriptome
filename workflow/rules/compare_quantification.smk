# unnormalised quantification correlation
rule quantification_getAverageCountsTables:
    input:
        ont_quant = [
            expand("data/quantification/teloprime/merged/teloprime_merged_{barcode}_quant.tsv",
                barcode=SAMPLE_INFO_ont["ont"]),
            expand("data/quantification/cdna/merged/cdna_merged_{barcode}_quant.tsv",
                barcode=SAMPLE_INFO_ont["cdna"]),
            expand("data/quantification/rna/merged/rna_merged_{barcode}_quant.tsv",
                barcode=["rt", "cool"])
            ],
        illumina_quant = expand("data/quantification/salmon/{sample}/quant.sf", sample=SAMPLES_ont),
        txdb = "data/annotation/annotation_txdb.sqlite"
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        gene = "data/comparisons/counts/comparisons_meanCounts_gene.tsv.gz",
        tx = "data/comparisons/counts/comparisons_meanCounts_tx.tsv.gz"
    script:
        "../scripts/comparisons_quantification_averageCountsTable.R"


rule quantification_correlation:
    input:
        biomaRt_gene = "data/annotation/biomaRt_gene.rds",
        biomaRt_tx = "data/annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        gene = "data/comparisons/counts/comparisons_meanCounts_gene.tsv.gz",
        tx = "data/comparisons/counts/comparisons_meanCounts_tx.tsv.gz"
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "res/comparisons/comparisons_quantification_correlation.html"
    script:
        "../scripts/comparisons_quantification_correlation.Rmd"


rule quantification_getCountsTables:
    input:
        ont_quant = [expand("data/quantification/teloprime/{flowcell}/teloprime_{flowcell}_{barcode}_quant.tsv",
                         flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["ont"]),
                     expand("data/quantification/cdna/{flowcell}/cdna_{flowcell}_{barcode}_quant.tsv",
                         flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["cdna"])],
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "data/comparisons/counts/comparisons_counts.tsv.gz"
    script:
        "../scripts/comparisons_quantification_countsTable.R"


rule quantification_correlationWithinSamples:
    input:
        "data/comparisons/counts/comparisons_counts.tsv.gz"
    output:
        "res/comparisons/comparisons_quantification_correlationWithinSamples.html"
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "../scripts/comparisons_quantification_correlationWithinSamples.Rmd"


# normalised quantification correlation
rule quantification_correlation_normalised:
    input:
        biomaRt_gene = "data/annotation/biomaRt_gene.rds",
        biomaRt_tx = "data/annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        gene = [expand("data/deseq/{dataset}/{dataset}_gene_cm_ntd.csv.gz",
                    dataset=["cdna", "rna", "teloprime"]),
                 "data/deseq/illumina/illumina_gene_cm_tpm.csv.gz"],
        tx = [expand("data/deseq/{dataset}/{dataset}_transcript_cm_ntd.csv.gz",
                  dataset=["cdna", "rna", "teloprime"]),
              "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz"]
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "res/comparisons/comparisons_quantification_correlation_normalised.html"
    script:
        "../scripts/comparisons_quantification_correlation_normalised.Rmd"


rule make_deseqDataSet_libraries:
    input:
        teloprime = expand("data/quantification/teloprime/{flowcell}/teloprime_{flowcell}_{barcode}_quant.tsv",
            flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["ont"]),
        cdna = expand("data/quantification/cdna/{flowcell}/cdna_{flowcell}_{barcode}_quant.tsv",
            flowcell=["flowcell1", "flowcell2"], barcode=SAMPLE_INFO_ont["cdna"]),
        sample_info = config["SAMPLE_INFO"]
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "data/comparisons/counts/counts_withinSamplesNorm.tsv.gz"
    script:
        "../scripts/comparisons_quantification_importONT_library.R"


rule quantification_correlationWithinSamples_normalised:
    input:
        "data/comparisons/counts/counts_withinSamplesNorm.tsv.gz"
    output:
        "res/comparisons/comparisons_quantification_correlationWithinSamples_normalised.html"
    conda:
        "../envs/r_4.1.2.yaml"
    script:
        "../scripts/comparisons_quantification_correlationWithinSamples_normalised.Rmd"


rule counts_pca:
    input:
        gene = [expand("data/deseq/{dataset}/{dataset}_gene_cm_ntd.csv.gz",
                    dataset=["cdna", "rna", "teloprime"]),
                 "data/deseq/illumina/illumina_gene_cm_tpm.csv.gz"],
        tx = [expand("data/deseq/{dataset}/{dataset}_transcript_cm_ntd.csv.gz",
                  dataset=["cdna", "rna", "teloprime"]),
              "data/deseq/illumina/illumina_transcript_cm_tpm.csv.gz"],
        sample_info = config["SAMPLE_INFO"]
    conda:
        "../envs/r_4.1.2.yaml"
    output:
        "res/comparisons/comparisons_countsPCA.html"
    script:
        "../scripts/comparisons_countsPCA.Rmd"
