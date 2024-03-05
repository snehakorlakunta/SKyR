#' @title threshold_cells
#'
#' @description Creates a logarithmic scatter plot of the given values
#'     with upper (gold) and lower (magenta) threshold lines on the y-axis
#'
#' @param title Name of dataset
#' @param values Values to be plotted on a log scale
#' @param min_per_cell Lower threshold value
#' @param max_per_cell Upper threshold value
#'
#' @return A log-scale scatterplot with upper/lower thresholds
#'
#' @export
threshold_cells <- function(title, values, min_per_cell, max_per_cell) {
  #counts per cell = nCount_RNA || genes per cell = nFeature_RNA
  plot(sort(values), xlab='cell', log='y', main=paste(title, "(ordered)"))
  abline(h=min_per_cell, col='magenta')  # lower threshold
  abline(h=max_per_cell, col='gold') # upper threshold
}
