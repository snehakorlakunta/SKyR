#' @title vln_merge_data
#'
#' @author Sneha Korlakunta
#'
#' @description Retrieves data from Seurat object for plotting multi-feature violin plots
#'              (only tested for 2 conditions).
#'
#' @param seu_obj A Seurat object
#' @param goi A list of genes of interest
#' @param group A column in the object's metadata to group the data by (legend)
#'
#' @return A table with Gene, Group, and Expression values.
#'
#' @export
#'
vln_merge_data <- function(seu_obj, goi, group = "CellType") {

  vln <- VlnPlot(seu_obj, goi, group.by = group)

  vln_data <- data.table::data.table()
  for (i in 1:length(goi)) {
    gene_data <- vln[[i]]$data
    gene_data$Gene <- goi[[i]]
    gene_data$Expression <- gene_data[[goi[[i]]]]
    gene_data[[goi[[i]]]] <- NULL

    vln_data <- rbind(gene_data, vln_data)
  }

  vln_data$Gene <- factor(vln_data$Gene, levels = goi)
  vln_data$ident <- factor(vln_data$ident)

  return(vln_data)
}
