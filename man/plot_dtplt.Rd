% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_dtplt.R
\name{plot_dtplt}
\alias{plot_dtplt}
\title{plot_dtplt}
\usage{
plot_dtplt(dtplt_data, split = FALSE)
}
\arguments{
\item{dtplt_data}{A data frame with the following columns: Gene, Grp, Split, Avg, Avg.SC, Pct, Pct.SC}

\item{split}{A logical value indicating the data is split by gene}
}
\value{
A ggplot object or a list of ggplot objects
}
\description{
Generates multi-feature plot(s).
}
