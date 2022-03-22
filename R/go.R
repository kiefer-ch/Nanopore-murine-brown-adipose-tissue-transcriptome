#' Do GO analysis using topGO
#'
#' `get.geneList()` makes a genelist for topgoi.
#'
#' @param significant Another list of ensembl gene identifiers.
#' @param reference A list of ensembl gene identifiers (without version).
#'
#' @return A named vector with 0 and 1.
#'
#' @export
get_geneList <- function(significant, reference) {

    checkmate::check_subset(significant, reference)

    gl <- factor(as.integer(reference %in% significant), levels = c(0, 1))
    names(gl) <- reference

    gl
}

#' `make.topGO()` makes a topGOdata object.
#'
#' @param geneList A geneList as made by get.geneList.
#' @param description Description.
#'
#' @return A topGOdata object.
#'
#' @export
make.topGO <- function(geneList, description) {
    new("topGOdata",
        description = description,
        ontology = "BP",
        allGenes = geneList,
        nodeSize = 10,
        annot = topGO::annFUN.org,
        ID = "ensembl",
        mapping = "org.Mm.eg")
}

#' `get.results()` gets results from a topGOdata object.
#'
#' @param topGO A topGOdata object.
#' @param method Algorithm to use within topGO::runTest.
#'
#' @return A tibble.
#'
#' @export
get.results <- function(topGO, method = "parentchild") {
    topGOres <- topGO::runTest(topGO, algorithm = method, statistic = "fisher")

    topGO::GenTable(topGO, p = topGOres, orderBy = "p", ranksOf = "p", topNodes = 100) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(p = dplyr::if_else(p == "< 1e-30", 10^-30, suppressWarnings(as.double(p)))) %>%
        dplyr::mutate(rank = dplyr::row_number())
}

#' `remove.version()` removes version from ensemble gene ids.
#'
#' @param x A  character vector with versions.
#'
#' @return A character vector without versions.
#'
#' @export
remove.version <- function(x) {
    unlist(lapply(stringr::str_split(x, "[.]"), "[[",1))
}

#' `plot.go.summary()` Plots a heatmap of p values keeping the n most
#'   significant go terms per group.
#'
#' @param res_go_list A list with individual results tibbles from get.results.
#' @param n Number of terms to keep per group.
#' @param n_char Max characters for terms. Otherwise cropped.
#' @param base_size Base size for text.
#' @param x_label_size Pch for x axis labels.
#'
#' @return A ggplot object.
#'
#' @export
plot.go.summary <- function(res_go_list, n = 5, n_char = 30, base_size = 8, x_label_size = 6) {
    # list of unique GO terms within the most n significant in at least one group
    go_to_keep <- res_go_list %>%
        purrr::map(dplyr::slice_min, rank, n = n) %>%
        purrr::map(dplyr::select, GO.ID, p) %>%
        dplyr::bind_rows(.id = "cluster") %>%
        dplyr::group_by(GO.ID) %>%
        dplyr::summarise(p = min(p, na.rm = TRUE)) %>%
        dplyr::pull(GO.ID)

    df_go <- res_go_list %>%
        purrr::map(mutate, rate = Significant/Annotated) %>%
        purrr::map(dplyr::select, GO.ID, p, rate) %>%
        dplyr::bind_rows(.id = "cluster") %>%
        dplyr::filter(GO.ID %in% go_to_keep) %>%
        tidyr::pivot_wider(
            names_from = cluster,
            values_from = c(p, rate),
            values_fill = list("p" = 1, "rate" = 0)) %>%
        dplyr::mutate(across(starts_with("p_"), log10))

    clust_go <- df_go %>%
        dplyr::select(matches("GO.ID|^p_")) %>%
        tibble::column_to_rownames("GO.ID") %>%
        data.matrix() %>%
        dist("euclidean") %>%
        hclust("complete")

    df_go %>%
        dplyr::mutate(GO.ID = factor(GO.ID, levels = df_go$GO.ID[clust_go$order])) %>%
        tidyr::pivot_longer(cols = matches("^p_|^rate_")) %>%
        tidyr::separate(name, c("type", "cluster"), sep = '_', extra = "merge") %>%
        tidyr::pivot_wider(names_from = type, values_from = value) %>%
        dplyr::mutate(p = 10^p) %>%
        mutate(rate = if_else(rate == 0, NA_real_, rate)) %>%
        ggplot2::ggplot(aes(cluster, GO.ID)) +
        ggplot2::geom_point(aes(colour = p, size = rate)) +
        ggplot2::scale_colour_viridis_c(direction = -1, trans = "log10") +
        ggthemes::theme_tufte(base_size = base_size, base_family = "Helvetica") +
        scale_y_discrete(name = NULL, breaks = df_go$GO.ID,
            labels = function(x) {Term(x) %>%
            dplyr::if_else(nchar(.) > n_char, paste0(substr(., 1, n_char - 3), "..."), .)}) +
        ggplot2::xlab("Cluster") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = x_label_size))
}

#' `plot.go()` removes version from ensemble gene ids.
#'
#' @param res_go A results tibble from get.results.
#' @param n Max number of elements to plot.
#' @param n_char Max number of chars for y axis labels.
#' @param base_size Base size for text.
#' @param alpha Alpha for geom_point.
#'
#' @return A ggplot object.
#'
#' @export
plot.go <- function(res_go, n = 15, n_char = 30, base_size = 8, alpha = .5,
    range = c(.1, 2.5)) {

    df_go <- res_go %>%
        dplyr::top_n(-n, rank) %>%
        dplyr::mutate(GeneRatio = Significant / Annotated) %>%
        dplyr::arrange(GeneRatio) %>%
        dplyr::mutate(GO.ID = factor(GO.ID, levels = GO.ID))

    df_go %>%
        ggplot2::ggplot(ggplot2::aes(GeneRatio, GO.ID)) +
        ggplot2::geom_point(ggplot2::aes(colour = p, size = Significant), alpha = alpha) +
        ggplot2::scale_color_viridis_c(direction = -1, trans = "log10") +
        ggplot2::scale_size_continuous(range = range) +
        ggthemes::theme_tufte(base_size = base_size, base_family = "Helvetica") +
        ggthemes::geom_rangeframe() +
        ggplot2::scale_y_discrete(name = NULL, breaks = df_go$GO.ID,
            labels = function(x) {Term(x) %>%
                dplyr::if_else(nchar(.) > n_char, paste0(substr(., 1, n_char - 3), "..."), .)})
}
