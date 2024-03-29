% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_duplicate_classification.R
\name{get_tandem_proximal}
\alias{get_tandem_proximal}
\title{Classify gene pairs derived from tandem and proximal duplications}
\usage{
get_tandem_proximal(pairs = NULL, annotation_granges = NULL, proximal_max = 10)
}
\arguments{
\item{pairs}{A 3-column data frame with columns \strong{dup1}, \strong{dup2},
and \strong{type} indicating duplicated gene 1, duplicated gene 2, and
the mode of duplication associated with the pair. This data frame
is returned by \code{get_segmental()}.}

\item{annotation_granges}{A processed GRanges object as in each element
of the list returned by \code{syntenet::process_input()}.}

\item{proximal_max}{Numeric scalar with the maximum distance (in number
of genes) between two genes to consider them as proximal duplicates.
Default: 10.}
}
\value{
A 3-column data frame with the variables:
\describe{
\item{dup1}{Character, duplicated gene 1.}
\item{dup2}{Character, duplicated gene 2.}
\item{type}{Factor of duplication types, with levels
"SD" (segmental duplication),
"TD" (tandem duplication),
"PD" (proximal duplication), and
"DD" (dispersed duplication).}
}
}
\description{
Classify gene pairs derived from tandem and proximal duplications
}
\examples{
data(yeast_annot)
data(yeast_seq)
data(fungi_kaks)
scerevisiae_kaks <- fungi_kaks$saccharomyces_cerevisiae

# Get processed annotation for S. cerevisiae
pdata <- annotation <- syntenet::process_input(yeast_seq, yeast_annot)
annot <- pdata$annotation[[1]]

# Get duplicated pairs
pairs <- scerevisiae_kaks[, c("dup1", "dup2", "type")]
pairs$dup1 <- paste0("Sce_", pairs$dup1)
pairs$dup2 <- paste0("Sce_", pairs$dup2)

# Get tandem and proximal duplicates
td_pd_pairs <- get_tandem_proximal(pairs, annot)

}
