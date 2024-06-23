#' @title remove_cells
#'
#' @description Allows user to manually remove cells in a Seurat object.
#'
#' @param seurat_object A Seurat object
#'
#' @return A Seurat object with selected cells removed.
#'
#' @export
#'
remove_cells <- function(seurat_object) {
  print("Cancel function at any point by selecting '0'.")

  select1 <- utils::menu(c("Whole object", "Subcluster(s)", "Whole object, but remove from subcluster(s)"), title = "Return whole object or a subcluster?")
  if (select1 == 0) {stop()}

  if (select1 %in% c(2, 3)) {
    subclusters <- levels(seurat_object$CellType)
    select_subcluster <- menu(c(subclusters, "Multiple"), title = "Which subcluster?")
    if (select_subcluster == 0) {stop()}
    else if (select_subcluster == length(subclusters) + 1) {
      print("Select multiple subtypes, separated by spaces. Input `enter` twice to complete operation.")
      all_subtypes <- scan()
      print("Selected subtypes:", subclusters[all_subtypes])
    }
    else {
      all_subtypes <- as.integer(select_subcluster)
    }

    if (select1 == 3) {whole_object <- seurat_object}
    seurat_object <- subset(seurat_object, CellType %in% subclusters[all_subtypes])
  }

  cells_to_remove <- select_cells(seurat_object)

  if (select1 == 3) {seurat_object_RM <- whole_object[,!colnames(whole_object) %in% cells_to_remove]}
  else {seurat_object_RM <- seurat_object[,!colnames(seurat_object) %in% cells_to_remove]}

  print(Seurat::DimPlot(seurat_object_RM))

  return (seurat_object_RM)
}
