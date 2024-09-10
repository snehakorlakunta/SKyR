#' @title get_dtplt_data
#'
#' @author Sneha Korlakunta
#'
#' @description Retrieves data from Seurat object for plotting multi-feature dot plots.
#'
#' @param seu_obj A Seurat object
#' @param goi A list of genes of interest
#' @param group A column in the object's metadata to group the data by (x-axis)
#' @param split A column in the object to split the data by (legend)
#' @param sep A separator for the group and split columns
#' @param scale_100 Scale the percentage data to 0-100
#'
#' @return A table with Gene, Group (Grp), Split (if applied), Average, scaled Average,
#'         Percentage, and scaled Percentage values.
#'
#' @export
#'
get_dtplt_data <- function(seu_obj, goi, group = "CellType", split = NULL, sep = "_", scale_100 = TRUE) {

  goi <- intersect(goi, rownames(seu_obj)) # ensures genes are available to plot in the object

  genes <- list()
  grp <- list()
  avg <- list()
  avg.sc <- list()
  pct <- list()
  pct.sc <- list()

  for (i in 1:length(goi)) {
    plt <- Seurat::DotPlot(seu_obj, goi[[i]], split.by = split, group.by = group, cols = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff'))

    plt.pct <- scale(x = log1p(plt$data$pct.exp))
    if (scale_100) {plt.pct <- MinMax(data = plt.pct, min = 0, max = 100)}
    else {plt.pct <- MinMax(data = plt.pct, min = min(plt.pct), max = max(plt.pct))}

    genes <- append(genes, plt$data$features.plot)
    grp <- append(grp, plt$data$id)
    avg <- append(avg, plt$data$avg.exp)
    avg.sc <- append(avg.sc, plt$data$avg.exp.scaled)
    pct <- append(pct, plt$data$pct.exp)
    pct.sc <- append(pct.sc, plt.pct)
  }

  if (!is.null(split)) {
    dt_data <- data.frame(
      Gene =       genes,
      Grp =        unlist(word(grp, 1, sep = sep)),
      Split =      unlist(word(grp, 2, sep = sep)),
      Avg =        unlist(avg, recursive = FALSE),
      Avg.SC =     unlist(avg.sc, recursive = FALSE),
      Pct =        unlist(pct, recursive = FALSE),
      Pct.SC =     unlist(pct.sc, recursive = FALSE)
    )
  } else {
    dt_data <- data.frame(
      Gene =       genes,
      Grp =        grp,
      Avg =        unlist(avg, recursive = FALSE),
      Avg.SC =     unlist(avg.sc, recursive = FALSE),
      Pct =        unlist(pct, recursive = FALSE),
      Pct.SC =     unlist(pct.sc, recursive = FALSE)
    )
  }

  return(dt_data)
}
