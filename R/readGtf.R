#' `read_gtf` reads gtf files.
#'
#' @param file A gtf file.
#'
#' @return A tibble.
#'
#' @export
#'
read_gtf <- function(file) {
    readr::read_delim(file = file,
        delim = "\t",
        comment = "#",
        na = c('.'),
        col_names = c("sequence", "source", "feature", "start", "end", "score",
            "strand", "phase", "attributes"),
        col_types = readr::cols(
            sequence = col_character(),
            source = col_character(),
            feature = col_character(),
            start = col_integer(),
            end = col_integer(),
            score = col_character(),
            strand = col_character(),
            phase = col_character(),
            attributes = col_character()),
        progress = FALSE) %>%
    dplyr::filter(strand %in% c('+', '-')) %>%
    dplyr::mutate(feature = as.factor(feature),
        strand = as.factor(strand))
}
