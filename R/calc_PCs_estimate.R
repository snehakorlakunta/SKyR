#' @title calc_PCs_estimate
#'
#' @description Calculates an estimate of how many PCs cover most of the variance in a dataset.
#'
#' @param seurat_object A seurat object
#' @param slot A dimensional reduction present in the seurat object
#'
#' @return minimum of co1 and co2 (and prints co1 and co2)
#'
#' @export
#'
calc_PCs_estimate <- function(seurat_object, slot = "pca") {
  #https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
  # Determine percent of variation associated with each PC
  pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]

  #co1
  #The first metric returns PC42 as the PC matching these requirements. Let’s check the second metric, which identifies the PC where the percent change in variation between consecutive PCs is less than 0.1%:

  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  # last point where change of % of variation is more than 0.1%.
  #co2
  #This second metric returns PC14. Usually, we would choose the minimum of these two metrics as the PCs covering the majority of the variation in the data.

  # Minimum of the two calculation
  pcs <- min(co1, co2)

  print(paste0("co1 = ", co1, ", co2 = ", co2))

  return(pcs)
}
