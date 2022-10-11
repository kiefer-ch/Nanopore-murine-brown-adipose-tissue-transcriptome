
source(".Rprofile")
suppressPackageStartupMessages({
    library("dplyr")
    library("readr")
    library("purrr")
})

################################################################################
#
# author: Christoph Kiefer
# email: christophak@bmb.sdu.dk
#
################################################################################

genes1 <- c("Gtf2b", "Hprt", "Cars2_long", "Cars2_short", "Scp2_5p", "Scp2_3p",
            "Ergic1_long", "Ergic1_short", "Smyd4_all", "Smyd4_long", "Adtrp_long",
            "Adtrp_all")
genes2 <- c("Pex6_long", "Pex6_all", "Mlxipl_long", "Mlxipl_short",
            "Dipk1b_short", "Dipk1b_long", "Acsl5_short", "Acsl5_long", "Gnas_long",
            "Gnas_short", "Aldoa_long", "Aldoa_short")
genes3 <- c("Adcy3_long", "Adcy3_short", "Lipe_long", "Lipe_short",
            "Ppargc1a_long", "Ppargc1a_short", "Pde4d_long", "Pde4d_short")
hkp <- c("Gtf2b", "Hprt")
notPass <- c()


sample_info <- read_csv(snakemake@input[["sample_info"]]) %>%
    mutate_at(vars(matches("condition")), as.factor) %>%
    select_at(vars(matches("sample_id|condition"))) %>%
    filter(!sample_id %in% c("rt", "cool"))

get.pipettingScheme <- function(genes, sampleInfo, rep_qpcr = 2, nrow = 16) {

    n_samples <- nrow(sampleInfo)
    n_genes <- length(genes)

    length_gene <- n_samples * rep_qpcr
    cols_per_gene <- ((length_gene - 1) %/% nrow) + 1

    samples_per_column <- rep(c(rep(nrow, (length_gene %/% nrow)),
                                length_gene %% nrow)[c(rep(nrow, (length_gene %/% nrow)),
                                                       length_gene %% nrow) != 0], n_genes)

    sampleInfo <- sampleInfo %>%
        tidyr::uncount(rep_qpcr) %>%
        replicate(n_genes, ., simplify = FALSE) %>%
        bind_rows()

    tibble(gene = rep(genes, each = length_gene),
           col = rep(1:(cols_per_gene * n_genes), samples_per_column),
           row = rep(1:length_gene %% nrow, n_genes)) %>%
        mutate(row = LETTERS[if_else(row == 0, nrow, row)]) %>%
        mutate(pos = paste0(row, col)) %>%
        bind_cols(., sampleInfo) %>%
        dplyr::select(-col, -row, gene)
}

import.LCcq <- function(pipettingScheme, file, decimal_mark = '.', maxCq = 40) {

    df <- read_tsv(file, locale = locale(decimal_mark = decimal_mark),
                   skip = 1,
                   col_types = cols(
                       Include = col_character(),
                       Color = col_integer(),
                       Pos = col_character(),
                       Name = col_character(),
                       Cp = col_double(),
                       Concentration = col_character(),
                       Standard = col_integer(),
                       Status = col_character())) %>%
        dplyr::select(cq = Cp, pos = Pos) %>%
        mutate(cq = if_else(cq < maxCq, cq, NA_real_))
    pipettingScheme %>%
        left_join(df, by = "pos")
}

df1 <- get.pipettingScheme(
    genes = genes1,
    sampleInfo = sample_info,
    nrow = 16) %>%
    import.LCcq(file = snakemake@input[["cq1"]],
                decimal_mark = '.') %>%
    mutate(cq = if_else(pos %in% notPass, NA_real_, cq))

df2 <- get.pipettingScheme(
    genes = genes2,
    sampleInfo = sample_info,
    nrow = 16) %>%
    import.LCcq(file = snakemake@input[["cq2"]],
                decimal_mark = '.') %>%
    mutate(cq = if_else(pos %in% notPass, NA_real_, cq))

df3 <- get.pipettingScheme(
    genes = genes3,
    sampleInfo = sample_info,
    nrow = 16) %>%
    import.LCcq(file = snakemake@input[["cq3"]],
                decimal_mark = '.') %>%
    mutate(cq = if_else(pos %in% notPass, NA_real_, cq)) %>%
    filter(sample_id != "190220_5_iBAT") # there was not enough cDNA left

df <- list(df1, df2, df3)

# export cq values
df %>%
    write_rds(snakemake@output$cq)


detect.outliers <- function(df, cutoff = 1) {
    out <- df %>%
        group_by_at(vars(matches("sample_id|gene"))) %>%
        mutate(error = max(cq, na.rm = TRUE) - min(cq, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(high_var = if_else(error > cutoff, TRUE, FALSE)) %>%
        dplyr::select(-error)

    if(any(out$high_var, na.rm = TRUE)) {
        cat(sprintf("%i samples set to NA due to high variance in qPCR
                replicates.\n",
                    sum(out$high_var, na.rm = TRUE) / 2))
    } else {
        cat("No samples with high variance detected.\n")
    }

    return(out)
}


df <- df %>%
    map(detect.outliers, cutoff = 1)


get.averageCq <- function(df) {
    df %>%
        mutate(cq = if_else(high_var, NA_real_, cq)) %>%
        group_by_at(vars(matches("sample_id|condition|gene"))) %>%
        summarise(cq = mean(cq, na.rm = TRUE)) %>%
        ungroup()
}

df <- df %>%
    bind_rows() %>%
    get.averageCq()


get.dCq <- function(df, hkp) {
    ref <- df %>%
        filter(gene %in% hkp) %>%
        group_by_at(vars(matches("sample_id"))) %>%
        summarise(housekeeper = mean(cq)) %>%
        ungroup()

    df %>%
        filter(!gene %in% hkp) %>%
        left_join(ref, by = colnames(df)[grep("sample_id", colnames(df))]) %>%
        mutate(dcq = cq - housekeeper) %>%
        dplyr::select(-cq, -housekeeper)
}

df <- df %>%
    get.dCq(hkp = hkp)


get.averageDCq <- function(df) {
    df %>%
        group_by_at(vars(matches("condition|gene|sample_id"))) %>%
        summarise(dcq = mean(dcq, na.rm = TRUE)) %>%
        ungroup()
}

df <- df %>%
    get.averageDCq()

df %>%
    saveRDS(snakemake@output$dcq)

# only those samples that are in the RNA seq
df %>%
    filter(sample_id %in% c("190220_2_iBAT", "190220_4_iBAT", "190220_9_iBAT",
            "190220_11_iBAT", "190220_14_iBAT", "190220_15_iBAT")) %>%
    dplyr::select(gene, condition_temp, dcq) %>%
    tidyr::separate(gene, c("gene", "isoform"), sep = '_') %>%
    tidyr::drop_na() %>%
    split(.$gene) %>%
    map(lm, formula = dcq ~ isoform * condition_temp) %>%
    map(summary.lm) %>%
    map(function(x) x$coefficients[16]) %>%
    unlist() %>%
    tibble::enframe(name = "mgi_symbol", value = "qpcr_p") %>%
    mutate(qpcr_padj = p.adjust(qpcr_p)) %>%
    mutate(sign = symnum(qpcr_padj,
                         cutpoints = c(0, .001, .01, .05, .1, 1),
                         symbols = c("***", "**", "*", ".", " "))) %>%
    write_csv(snakemake@output$stats)
