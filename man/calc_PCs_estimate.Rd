% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_PCs_estimate.R
\name{calc_PCs_estimate}
\alias{calc_PCs_estimate}
\title{calc_PCs_estimate}
\usage{
calc_PCs_estimate(seurat_object, slot = "pca")
}
\arguments{
\item{seurat_object}{A Seurat object}

\item{slot}{A dimensional reduction present in the seurat object}
}
\value{
minimum of co1 and co2 (and prints co1 and co2)
}
\description{
Calculates an estimate of how many PCs cover most of the variance in a dataset.
}
