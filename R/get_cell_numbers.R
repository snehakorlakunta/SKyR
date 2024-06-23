#' @title get_cell_numbers
#'
#' @description Generates a table of cell numbers among given columns.
#'
#' @param seurat_object A Seurat object
#' @param columns a list of one or more metadata columns (default is idents)
#' @param order the metadata column by which to order the columns, must be numerical
#' @param decr_order A bool defining whether sorted order is decreasing (default is TRUE)
#'
#' @return a data table containing the counts of cells per given column
#'
#' @export
#'
get_cell_numbers <- function(seurat_object, columns = "idents", order = NULL, decr_order = TRUE) {
  if (columns == "idents") {
    seurat_object$col_of_interest <- SeuratObject::Idents(seurat_object)
    columns <- c("col_of_interest")
  }

  metadata <- seurat_object@meta.data %>% data.table::as.data.table
  counts <- metadata[, data.table::.N, by = columns]

  ## with additional casting after the counting
  #md[, .N, by = c("Sample", "seurat_clusters")] %>% dcast(., Sample ~ seurat_clusters, value.var = "N")

  if (!(is.na(order))) {
    counts <- counts[order(counts[[order]], decreasing = decr_order), ]
  }

}
