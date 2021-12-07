
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

log_info("Reading ONT transcript quantification...")
df <- snakemake@input$ont_quant %>%
    set_names(basename(.) %>%
                  sub("_quant.tsv", "", .)) %>%
    map(read_tsv, col_types = "ci") %>%
    bind_rows(.id = "id") %>%
    tidyr::separate("id", c("library", "flowcell", "barcode"))


log_info("Writing to disc...")
df %>%
    write_tsv(snakemake@output[[1]])


log_success("Done.")
