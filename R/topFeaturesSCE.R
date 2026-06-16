#' @importFrom SummarizedExperiment rowData
#' @importFrom MultiAssayExperiment metadata

topFeaturesSCE <- function(sce, sort.by = NULL, n = 10) {
    rd <- as.data.frame(rowData(sce))
    comparisons <- metadata(sce)$comparisons

    ## Detect unadjusted output by presence of _GAM/_msqrob2/_aggregated suffixes
    is_unadjusted <- any(grepl("_GAM$|_msqrob2$|_aggregated$", colnames(rd)))

    if (is.null(sort.by)) {
        if (is_unadjusted) {
            stop(
                "Output from adjust = 'none' detected: please specify ",
                "'sort.by' explicitly (e.g. sort.by = 'omnibus_aggregated')."
            )
        }
        sort.by <- if (length(comparisons) == 1L) comparisons[1L] else "omnibus"
    }

    if (!sort.by %in% colnames(rd)) {
        stop(
            "'sort.by' column '", sort.by, "' not found in rowData. ",
            "Available p-value columns: ",
            paste(c(comparisons, if (!is.null(comparisons)) "omnibus"), collapse = ", ")
        )
    }

    rd_sorted <- rd[order(rd[[sort.by]], na.last = TRUE), , drop = FALSE]
    head(rd_sorted, n = n)
}
