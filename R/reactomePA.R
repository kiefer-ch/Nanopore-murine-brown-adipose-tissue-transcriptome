#' Do Reactome analysis using ReactomePA
#'
#' `plot.reactome.summary()` Plots a heatmap of p values keeping the n most
#'   significant reactome pathways per group.
#'
#' @param res_list A list with individual results tibbles from get.results.
#' @param n Number of terms to keep per group.
#' @param n_char Max characters for terms. Otherwise cropped.
#' @param base_size Base size for text.
#' @param x_label_size Pch for x axis labels.
#'
#' @return A ggplot object.
#'
#' @export
plot.reactome.summary <- function(res_list, n = 5, n_char = 30, base_size = 8, x_label_size = 6) {
    # list of unique reactome pathways within the most n significant in at least one group
    to_keep <- res_list %>%
        purrr::map(`@`, "result") %>%
        purrr::map(tibble::as_tibble) %>%
        purrr::map(dplyr::select, Description, qvalue) %>%
        purrr::map(dplyr::slice_min, qvalue, n = n, with_ties = FALSE) %>%
        dplyr::bind_rows(.id = "cluster") %>%
        dplyr::filter(qvalue < .05) %>%
        dplyr::pull(Description)

    df_reactome <- res_list %>%
        purrr::map(`@`, "result") %>%
        purrr::map(dplyr::as_tibble) %>%
        purrr::map(dplyr::select, Description, qvalue)

    # remove empty lists
    keep <- df_reactome %>%
        purrr::map(pull, qvalue) %>%
        purrr::map(~ . < .05) %>%
        purrr::map(any) %>%
        unlist()

    df_reactome <- df_reactome[keep]

    # stop if list is empty
    if(purrr::is_empty(df_reactome)) {
        message("No significant reactome terms.")
    } else {
        df_reactome <- df_reactome %>%
            dplyr::bind_rows(.id = "cluster") %>%
            dplyr::filter(Description %in% to_keep) %>%
            dplyr::mutate(
                Description = dplyr::if_else(nchar(Description) > 40,
                    paste0(strtrim(Description, 37), "..."), Description)) %>%
            tidyr::spread(cluster, qvalue) %>%
            dplyr::mutate_if(is.numeric, function(x) tidyr::replace_na(x, 1)) %>%
            dplyr::mutate_if(is.numeric, log10)

        clust_reactome <- df_reactome %>%
            tibble::column_to_rownames("Description") %>%
            data.matrix() %>%
            dist("euclidean") %>%
            hclust("complete")

        df_reactome %>%
            mutate(Description = factor(Description, levels = df_reactome$Description[clust_reactome$order])) %>%
            tidyr::gather("cluster", "qvalue", -Description) %>%
            mutate(qvalue = 10^qvalue) %>%
            ggplot(aes(cluster, Description)) +
            geom_tile(aes(fill = qvalue)) +
            scale_fill_viridis_c(direction = -1, trans = "log10") +
            theme_tufte(base_size = 8, base_family = "Helvetica") +
            ylab(NULL) +
            xlab("Cluster") +
            theme(axis.text.y = element_text(size = 8))
    }
}
