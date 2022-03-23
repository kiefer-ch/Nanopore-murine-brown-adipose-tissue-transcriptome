#' `subset.txdb` subsets txdb objects to all information of a single gene.
#'
#' @param txdb txdb object.
#' @param txdb_dump dumped txdb object.
#' @param gene_id Gene is to subset to.
#'
#' @return A txdb object.
#'
#' @export
#'
subset.txdb <- function(txdb, txdb_dump, gene_id) {

    all_tx <- AnnotationDbi::select(
            txdb,
            keys = gene_id,
            columns = c("GENEID", "TXNAME"),
            keytype = "GENEID") %>%
        dplyr::pull(TXNAME)

    all_tx_ids <- txdb_dump[[1]] %>%
        tibble::as_tibble() %>%
        dplyr::filter(tx_name %in% all_tx) %>%
        dplyr::pull(tx_id)

    new_txdb_dump <- vector("list", 4)
    new_txdb_dump[1:3] <- txdb_dump[1:3] %>%
        purrr::map(dplyr::filter, tx_id %in% all_tx_ids)
    new_txdb_dump[4] <- txdb_dump[4]

    do.call(makeTxDb, new_txdb_dump)
}
