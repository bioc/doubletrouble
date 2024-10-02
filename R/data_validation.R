
#' Check if gene names in set 1 are present in set 2
#' 
#' @param ref_ids Character vector of reference gene set.
#' @param test_ids Character vector of test gene set.
#' @param setnames Character vector of length with set names. 
#' Default: \code{c("gene pairs", "CDS")}
#' 
#' @return TRUE if names match, otherwise an error is shown.
#' @importFrom utils head
#' @details 
#' This internal function can be used, for instance, to check if CDS names
#' match gene IDs in the gene pair list.
#' @noRd
check_geneid_match <- function(
        ref_ids, test_ids, setnames = c("gene pairs", "CDS")
) {
    
    mismatch_ids <- ref_ids[!ref_ids %in% test_ids]
    mismatch_perc <- length(mismatch_ids) / length(ref_ids)
    mismatch_perc <- round(mismatch_perc * 100, 2)
    
    if(mismatch_perc >0) {
        stop(
            mismatch_perc, "%", " (N=", length(mismatch_ids), ") of the IDs in ", setnames[1], 
            " were not found in ", setnames[2], ".\n", 
            "All gene IDs in ", setnames[1], " must be in ", setnames[2], 
            ". Did you check if gene IDs match?",
            "\n\nHere are some examples of nonmatching IDs (from ", setnames[1], ") :\n",
            paste0(head(mismatch_ids, n = 5), collapse = "\n"),
            "\n\nAnd here are some examples of IDs in ", setnames[2], ":\n",
            paste0(head(test_ids, n = 5), collapse = "\n")
        )
    }
    
    return(TRUE)
}