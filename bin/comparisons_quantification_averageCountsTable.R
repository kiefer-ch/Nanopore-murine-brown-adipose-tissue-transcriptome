
source(".Rprofile")
suppressPackageStartupMessages({
    library("logger")
    library("dplyr")
    library("readr")
    library("purrr")
    library("tximport")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

log_info("Reading transcript to gene mapping...")
tx2g <- AnnotationDbi::loadDb(snakemake@input$"txdb")  %>%
    AnnotationDbi::select(.,
                          keys =  AnnotationDbi::keys(., keytype = "GENEID"),
                          keytype = "GENEID",
                          columns = "TXNAME") %>%
    dplyr::select(TXNAME, GENEID)


log_info("Reading ONT transcript quantification...")
df_ont_tx <- snakemake@input$ont_quant %>%
    set_names(basename(.) %>%
                  sub("_merged.+", "", .)) %>%
    map(read_tsv, col_types = "ci") %>%
    bind_rows(.id = "library") %>%
    group_by(library, Name) %>%
    summarise(NumReads = expm1(mean(log1p(NumReads)))) %>%
    tidyr::pivot_wider(names_from = library,
        values_from = NumReads,
        values_fill = 0L)


log_info("Reading illumina transcript quantification...")
df_illuma_tx <- tximport(files = snakemake@input$illumina_quant,
        type = "salmon",
        tx2gene = tx2g,
        txOut = TRUE,
        dropInfReps = TRUE)


log_info("Aggregating to gene level...")
df_illuma_gene <- summarizeToGene(
        df_illuma_tx,
        tx2g)

df_ont_gene <- df_ont_tx %>%
    left_join(tx2g, by = c("Name" = "TXNAME")) %>%
    group_by(GENEID) %>%
    summarise_if(is.numeric, sum)


log_info("Combining ONT and Illumina...")
df_tx <- df_illuma_tx$abundance %>%
    as_tibble(rownames = "ensembl_transcript_id_version") %>%
    mutate_at(vars(starts_with("V")), log1p) %>%
    mutate(illumina = expm1(rowMeans(dplyr::select(., starts_with("V"))))) %>%
    dplyr::select(ensembl_transcript_id_version, illumina) %>%
    full_join(df_ont_tx, by = c("ensembl_transcript_id_version" = "Name")) %>%
    mutate_if(is.numeric, tidyr::replace_na, 0) %>%
    filter(rowSums(dplyr::select_if(., is.numeric))  > 0)

df_gene <- df_illuma_gene$abundance %>%
    as_tibble(rownames = "ensembl_gene_id_version") %>%
    mutate_at(vars(starts_with("V")), log1p) %>%
    mutate(illumina = expm1(rowMeans(dplyr::select(., starts_with("V"))))) %>%
    dplyr::select(ensembl_gene_id_version, illumina) %>%
    full_join(df_ont_gene, by = c("ensembl_gene_id_version" = "GENEID")) %>%
    mutate_if(is.numeric, tidyr::replace_na, 0) %>%
    filter(rowSums(dplyr::select_if(., is.numeric))  > 0)


log_info("Combining ONT and Illumina...")
df_tx %>%
    write_tsv(snakemake@output$tx)

df_gene %>%
    write_tsv(snakemake@output$gene)


log_success("Done.")
