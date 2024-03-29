% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_ks_distro}
\alias{plot_ks_distro}
\title{Plot distribution of synonymous substitution rates (Ks)}
\usage{
plot_ks_distro(
  ks_df,
  min_ks = 0.01,
  max_ks = 2,
  bytype = FALSE,
  type_levels = NULL,
  plot_type = "histogram",
  binwidth = 0.03
)
}
\arguments{
\item{ks_df}{A data frame with Ks values for each gene pair
as returned by \code{pairs2kaks()}.}

\item{min_ks}{Numeric indicating the minimum Ks value to keep.
Default: 0.01.}

\item{max_ks}{Numeric indicating the maximum Ks value to keep.
Default: 2.}

\item{bytype}{Logical indicating whether or not to plot the distribution
by type of duplication (requires a column named \code{type}).}

\item{type_levels}{(Only valid if \strong{bytype} is not NULL) Character
indicating which levels of the variable specified in
parameter \strong{group_by} should be kept. By default, all levels are kept.}

\item{plot_type}{Character indicating the type of plot to create.
If \strong{bytype = TRUE}, possible types are "histogram" or "violin".
If \strong{bytype = FALSE}, possible types are "histogram", "density",
or "density_histogram". Default: "histogram".}

\item{binwidth}{(Only valid if \strong{plot_type = "histogram"})
Numeric indicating the bin width. Default: 0.03.}
}
\value{
A ggplot object.
}
\description{
Plot distribution of synonymous substitution rates (Ks)
}
\examples{
data(fungi_kaks)
ks_df <- fungi_kaks$saccharomyces_cerevisiae

# Plot distro
plot_ks_distro(ks_df, bytype = TRUE)
}
