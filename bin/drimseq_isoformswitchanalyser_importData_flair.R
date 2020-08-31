#!/usr/bin/Rscript --no-restore --no-environ --no-save

source(".Rprofile")
library("readr")
library("dplyr")
library("IsoformSwitchAnalyzeR")

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

# prepare sample_info
sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    filter(!is.na(ont)) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    mutate(path = file.path("salmon", illumina, "quant.sf"))

design_matrix <- sample_info %>%
    dplyr::select(sampleID = "sample_id",
        condition = condition_temp) %>%
    as.data.frame()

# import counts
counts <- read_tsv(snakemake@input[["counts"]]) %>%
    rename_at(vars(matches("condition")), function(x) strsplit(x, "_cond") %>% purrr::map(`[`, 1)) %>%
    tidyr::separate(ids, "isoform_id", sep = '_', extra = "drop") %>%
    as.data.frame()

# import drimseq results
test <- read_csv(snakemake@input[["test"]])

# crate switchlist
switchList <- importRdata(
    isoformCountMatrix   = counts,
    designMatrix         = design_matrix,
    isoformExonAnnoation = snakemake@input[["gtf"]],
    isoformNtFasta       = snakemake@input[["transcripts"]])

# add adjusted pvalues from drimseq
gene_lookup <- test$gene
names(gene_lookup) <- test$ensembl_gene_id_version

tx_lookup <- test$transcript
names(tx_lookup) <- test$ensembl_transcript_id_version

switchList$isoformFeatures$gene_switch_q_value <- gene_lookup[switchList$isoformFeatures$gene_id]
switchList$isoformFeatures$isoform_switch_q_value <- tx_lookup[switchList$isoformFeatures$isoform_id]

# annotate ORFs
switchList <- analyzeORF(switchList)

# add aa sequence
switchList <- extractSequence(switchList,
    extractNTseq = FALSE,
    writeToFile = TRUE,
    pathToOutput = snakemake@params[["aa_path"]],
    outputPrefix = snakemake@params[["aa_prefix"]])

# predict alternative splicing and intron retention
switchList <- analyzeAlternativeSplicing(switchList)
switchList <- analyzeIntronRetention(switchList)

# orf seq similarity
switchList <- analyzeSwitchConsequences(switchList,
    c("tss", "tts", "intron_retention", "ORF_seq_similarity", "NMD_status"))

saveRDS(switchList, snakemake@output[[1]])
