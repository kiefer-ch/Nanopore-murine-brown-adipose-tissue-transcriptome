# comparisons
rule quantification_averageCounts_tables:
    input:
        tx_counts = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_ntd.csv.gz",
                           dataset=["cdna", "teloprime", "rna"]),
        tx_tpm = "res/deseq/illumina/txlevel/illumina_txlevel_cm_ntd.csv.gz",
        gene_counts = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_cm_ntd.csv.gz",
                             dataset=["cdna", "teloprime", "rna"]),
        gene_tpm = "res/deseq/illumina/genelevel/illumina_genelevel_cm_ntd.csv.gz",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv"
    output:
        gene_avg = "res/comparisons/comparisons_meanCounts_gene.csv.gz",
        tx_avg = "res/comparisons/comparisons_meanCounts_tx.csv.gz",
        gene = "res/comparisons/comparisons_counts_gene.csv.gz",
        tx = "res/comparisons/comparisons_counts_tx.csv.gz"
    script:
        "comparisons_quantification_averageCountsTable.R"


rule quantification_correlation:
    input:
        gene_avg = "res/comparisons/comparisons_meanCounts_gene.csv.gz",
        tx_avg = "res/comparisons/comparisons_meanCounts_tx.csv.gz",
        gene = "res/comparisons/comparisons_counts_gene.csv.gz",
        tx = "res/comparisons/comparisons_counts_tx.csv.gz"
    output:
        "res/comparisons/comparisons_quantification_correlation.html"
    script:
        "comparisons_quantification_correlation.Rmd"


rule counts_pca:
    input:
        tx_counts = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
                           dataset=["cdna", "teloprime", "illumina", "rna"]),
        gene_counts = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_cm_cts.csv.gz",
                             dataset=["cdna", "teloprime", "illumina", "rna"]),
        sample_info = "sample_info/sampleInfo.csv"
    output:
        "res/comparisons/comparisons_countsPCA.html"
    script:
        "comparisons_countsPCA.Rmd"


rule feature_detection:
    input:
        tx_counts = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
                           dataset=["cdna", "teloprime", "illumina", "rna"]),
        gene_counts = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_cm_cts.csv.gz",
                             dataset=["cdna", "teloprime", "illumina", "rna"]),
        gene_tpm = "res/deseq/illumina/genelevel/illumina_genelevel_cm_ntd.csv.gz",
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        biotype_groups = "data/biotype_groups.csv",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        fig_folder = "res/fig/feature_detection"
    output:
        "res/comparisons/comparisons_feature_detection.html"
    script:
        "comparisons_featureDetection.Rmd"


rule bam_getAlignedLength:
    input:
        "bam/{dataset}/{file}_{type}.bam",
        "bam/{dataset}/{file}_{type}.bam.bai"
    wildcard_constraints:
        type = "genome|transcriptome",
        datatset = "teloprime|cdna"
    output:
        "res/comparisons/countReads/{dataset}_{file}_bam_{type}.rds"
    script:
        "comparisons_bam_getAlignedLength.R"


rule bam_getAlignedLength_rna:
    input:
        "bam/rna/{type}_{temp}_q7_sort.bam",
        "bam/rna/{type}_{temp}_q7_sort.bam.bai"
    wildcard_constraints:
        type = "genome|transcriptome"
    output:
        "res/comparisons/countReads/rna_{temp}_bam_{type}.rds"
    script:
        "comparisons_bam_getAlignedLength.R"


def get_bamnames(wildcards):
    files = list()
    if wildcards.dataset == "teloprime":
        for barcode in SAMPLE_INFO_ont["ont"]:
            for flowcell in ["X1_flowcell", "X3_flowcell"]:
                shortFlowcell = flowcell[0:2]
                for type in ["transcriptome", "genome"]:
                    filename = "bam/{}/{}/20191107_{}_{}_{}_q7_sort.bam".format(
                        wildcards.dataset, flowcell, shortFlowcell, type, barcode)
                    files.append(filename)
    elif wildcards.dataset == "cdna":
        for barcode in SAMPLE_INFO_ont["cdna"]:
            for pool in ["pool1", "pool2"]:
                for type in ["transcriptome", "genome"]:
                    filename = "bam/{}/{}/20200108_{}_{}_{}_q7_sort.bam".format(
                        wildcards.dataset, pool, pool, type, barcode)
                    files.append(filename)
    elif wildcards.dataset == "rna":
        for temperature in ["cool", "rt"]:
            for type in ["transcriptome", "genome"]:
                filename = "bam/rna/{}_{}_q7_sort.bam".format(
                    type, temperature)
                files.append(filename)
    return files


rule bam_getFlagStats_ont:
    input:
        get_bamnames
    output:
        "res/comparisons/countReads/{dataset}_flagstats.csv"
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    script:
        "comparisons_bam_countFlags.py"


def get_input_bam_coverage(wildcards):
    if wildcards.dataset in ["teloprime", "cdna"]:
        file_name = "bam/{}/{}_transcriptome.bam".format(
            wildcards.dataset, wildcards.barcode)
    elif wildcards.dataset == "rna":
        file_name = "bam/rna/transcriptome_{}_q7_sort.bam".format(wildcards.barcode)
    return [file_name, file_name + ".bai"]


rule bam_getCoverage:
    input:
        get_input_bam_coverage
    wildcard_constraints:
        dataset = "teloprime|cdna|rna"
    output:
        "res/comparisons/coverage/{dataset}/{barcode}.rds"
    script:
        "comparisons_bam_getCoverage.R"


rule coverage:
    input:
        geneBodyCoverage_teloprime = expand("res/comparisons/geneBody_coverage/teloprime/{barcode}.geneBodyCoverage.txt",
                                            barcode=SAMPLE_INFO_ont["ont"]),
        geneBodyCoverage_cdna = expand("res/comparisons/geneBody_coverage/cdna/{barcode}.geneBodyCoverage.txt",
                                       barcode=SAMPLE_INFO_ont["cdna"]),
        geneBodyCoverage_illumina = expand("qc/RSeQC/geneBody_coverage/{sample}.geneBodyCoverage.txt",
                                           sample=SAMPLES_ont),
        geneBodyCoverage_rna = expand("res/comparisons/geneBody_coverage/rna/{barcode}.geneBodyCoverage.txt",
                                       barcode=["rt", "cool"]),
        coverage_teloprime = expand("res/comparisons/coverage/teloprime/{barcode}.rds",
                                    barcode=SAMPLE_INFO_ont["ont"]),
        coverage_cdna = expand("res/comparisons/coverage/cdna/{barcode}.rds",
                               barcode=SAMPLE_INFO_ont["cdna"]),
        coverage_rna = expand("res/comparisons/coverage/rna/{barcode}.rds",
                               barcode=["rt", "cool"]),
        sample_info = "sample_info/sampleInfo.csv",
        biomaRt_tx = "annotation/biomaRt_tx.rds"
    params:
        fig_folder = "res/fig/coverage"
    output:
        "res/comparisons/comparisons_coverage.html"
    script:
        "comparisons_coverage.Rmd"


rule compare_dtu_1:
    input:
        drimseq = expand("res/drimseq/{dataset}/{dataset}_drimSeqStageR.csv",
                         dataset=["illumina", "cdna", "teloprime", "cdna_flair", "teloprime_flair"]),
        dexseq = expand("res/dexseq/{dataset}/{dataset}_dexseq_results.csv.gz",
                        dataset=["illumina", "cdna", "teloprime"]),
        biomaRt_gene = "annotation/biomaRt_gene.rds",
        counts = "res/comparisons/comparisons_meanCounts_gene.csv.gz"
    output:
        signif = "res/comparisons/comparisons_dtu_significant.csv",
        all = "res/comparisons/comparisons_dtu_all.rds"
    script:
        "comparisons_dtu.R"

rule compare_dtu_2:
    input:
        signif = "res/comparisons/comparisons_dtu_significant.csv",
        all = "res/comparisons/comparisons_dtu_all.rds"
    output:
        "res/comparisons/comparisons_dtu.html"
    script:
        "comparisons_dtu.Rmd"


def get_fastqnames(wildcards):
    files = list()
    if wildcards.dataset == "teloprime":
        for barcode in SAMPLE_INFO_ont["ont"]:
            for flowcell in ["X1_flowcell", "X3_flowcell"]:
                filename = "fastq/{}/{}/{}_q7.fastq.gz".format(
                    wildcards.dataset, flowcell, barcode)
                files.append(filename)
    elif wildcards.dataset == "cdna":
        for barcode in SAMPLE_INFO_ont["cdna"]:
            for pool in ["pool1", "pool2"]:
                filename = "fastq/{}/{}/{}_q7.fastq.gz".format(
                    wildcards.dataset, pool, barcode)
                files.append(filename)
    elif wildcards.dataset == "rna":
        for temperature in ["rt", "cool"]:
            filename = "fastq/{}/{}_q7.fastq.gz".format(
                wildcards.dataset, temperature)
            files.append(filename)
    return files


rule fastq_readLengthHistogram:
    input:
        get_fastqnames
    output:
        "res/comparisons/readLengthDistribution/{dataset}_fastqReadLengths.csv"
    script:
        "comparisons_fastq_readLengthHistogram.py"


rule read_lengths_fastq:
    input:
        readLengths = expand("res/comparisons/readLengthDistribution/{dataset}_fastqReadLengths.csv",
                             dataset=["teloprime", "cdna", "rna"]),
        avg_counts = "res/comparisons/comparisons_meanCounts_tx.csv.gz",
        sample_info = "sample_info/sampleInfo.csv"
    params:
        fig_folder = "res/fig/read_lengths"
    output:
        "res/comparisons/comparisons_readLengths_fastq.html"
    script:
        "comparisons_read_lengths_fastq.Rmd"


rule read_lengths_bam:
    input:
        teloprime_bam_tx = expand("res/comparisons/countReads/teloprime_{barcode}_bam_transcriptome.rds",
                                  barcode=SAMPLE_INFO_ont["ont"]),
        cdna_bam_tx = expand("res/comparisons/countReads/cdna_{barcode}_bam_transcriptome.rds",
                             barcode=SAMPLE_INFO_ont["cdna"]),
        biomaRt_tx = "annotation/biomaRt_tx.rds",
        sample_info = "sample_info/sampleInfo.csv",
        flagstats = expand("res/comparisons/countReads/{dataset}_flagstats.csv",
                           dataset=["cdna", "teloprime"])
    params:
        fig_folder = "res/fig/read_lengths"
    output:
        "res/comparisons/comparisons_readLengths_bam.html"
    script:
        "comparisons_read_lengths_bam.Rmd"


rule compare_dge:
    input:
        dte = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_apeglm_results.csv.gz",
                     dataset=["illumina", "teloprime", "cdna"]),
        dge = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_apeglm_results.csv.gz",
                     dataset=["illumina", "teloprime", "cdna"]),
        gene_counts = expand("res/deseq/{dataset}/genelevel/{dataset}_genelevel_cm_cts.csv.gz",
                             dataset=["cdna", "teloprime", "illumina"]),
        tx_counts = expand("res/deseq/{dataset}/txlevel/{dataset}_txlevel_cm_cts.csv.gz",
                           dataset=["cdna", "teloprime", "illumina"]),
        grouped_biotypes = "data/biotype_groups.csv"
    output:
        "res/comparisons/comparisons_dgeDte.html"
    script:
        "comparisons_dgeDte.Rmd"


rule compare_reannotation:
    input:
        gffcompare = ["stringtie/illumina/all_gffcompare.stringtie_illumina_merged.gtf.tmap",
                      expand("flair/{dataset}/all_gffcompare.flair.collapse.isoforms.gtf.tmap",
                             dataset=["cdna", "teloprime"])],
        sqanti = ["stringtie/illumina/sqanti/stringtie_illumina_merged_noUnknownStrand_classification.txt",
                  expand("flair/{dataset}/sqanti/flair.collapse.isoforms_classification.txt",
                         dataset=["cdna", "teloprime"])]
    output:
        "res/comparisons/comparisons_reannotation.html"
    script:
        "comparisons_reannotation.Rmd"


rule render_GOcomp:
    input:
        go = expand("res/deseq/{dataset}/{level}/{dataset}_{level}_apeglm_topgo.rds",
            dataset=["illumina", "cdna", "teloprime"],
            level=["txlevel", "genelevel"])
    params:
        cutoff = .01
    output:
        "res/comparisons/comparisons_go.html"
    script:
        "comparisons_go.Rmd"
