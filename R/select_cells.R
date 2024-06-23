#' @title select_cells
#'
#' @description Allows user to manually select cells in a Seurat object.
#'
#' @param seurat_object A Seurat object
#' @param selected_cells An internal parameter
#'
#' @return A list of cells selected from the input Seurat object.
#'
#' @export
#'
select_cells <- function(seurat_object, selected_cells = list()) {
  result <- Seurat::CellSelector(Seurat::DimPlot(seurat_object))   #https://www.biostars.org/p/9582411/#9582419
  if (length(selected_cells) != 0) {result <- c(result, selected_cells)}

  select2 <- utils::menu(c("Looks good!", "Select more cells", "Redo selection"), title = "How would you like to proceed?")
  if (select2 == 0) {stop()}

  if (select2 == 1) {
    print(Seurat::DimPlot(seurat_object, cells.highlight = result, cols.highlight = "darkred", cols= "grey"))
    return (result)
  }
  else if (select2 == 2) {
    select_cells(seurat_object, selected_cells = result)
  }
  else if (select2 == 3) {
    select_cells(seurat_object)
  }
}
