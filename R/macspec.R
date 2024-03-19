#' @title macspec
#'
#' @description Wrapper function to calculate MPI and AMDI of macrophage samples
#'     using MacSpectrum's algorithm.
#'
#' @param main a Seurat object containing myeloid cells
#' @param myeloid_sub subset of 'main' containing only myeloid/macrophage cells
#' @param myeloid_compare column of interest in Seurat object
#' @param select_hu_mo a string indicating if the samples are of human ("hum") or mouse ("mou"). default "mou"
#' @param v3_assay a string indicating which assay to retrieve from a Seurat v3/v4 object. default "SCT"
#' @param v5_layer a string indicating which layer to retrieve from a Seurat v5 object. default "data"
#'
#' @return Seurat object 'main' with a new metadata column called 'CellType_macspec'
#'
#' @export
macspec <- function(main, myeloid_sub, myeloid_compare, select_hu_mo="mou", v3_assay = "SCT", v5_layer = "data") {
  if (substr(SeuratObject::Version(myeloid_sub), 1, 1) == 5) {
    expression_matrix <- SeuratObject::LayerData(object = myeloid_sub, layer = v5_layer)
    expression_matrix <- as.matrix(expression_matrix)
    myeloid_mat <- as(expression_matrix, "sparseMatrix")
  } else {
    myeloid_mat <- myeloid_sub@assays[[v3_assay]]$data  # Seurat v3/4
  }

  rm(expression_matrix)

  myeloid_mat <- as.matrix(myeloid_mat)
  myeloid_df <- as.data.frame((myeloid_mat))
  myeloid_df <- myeloid_df %>% tibble::rownames_to_column(var="geneid")

  if (select_hu_mo == "mou") {
    myeloid_df[, 1] <- AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
                              keys = myeloid_df[, 1],
                              column = "ENSEMBL",
                              keytype = "SYMBOL",
                              multiVals = "first")
  } else {
    myeloid_df[, 1] <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                              keys = myeloid_df[, 1],
                              column = "ENSEMBL",
                              keytype = "SYMBOL",
                              multiVals = "first")
  }

  print(sum(rowSums(is.na(myeloid_df))))

  myeloid_df_na <- na.omit(myeloid_df)

  myeloid_class <- uconn_macspectrum(myeloid_df_na, myeloid_compare, select_hu_mo = select_hu_mo)

  matching_indices <- match(rownames(myeloid_sub@meta.data), rownames(myeloid_class))
  myeloid_sub@meta.data$MPI <- myeloid_class$MPI[matching_indices]
  myeloid_sub@meta.data$AMDI <- myeloid_class$AMDI[matching_indices]

  myeloid_sub@meta.data <- myeloid_sub@meta.data %>%
    dplyr::mutate(macrophage = dplyr::case_when(
      AMDI < 0 & MPI < 0 ~ "M0",
      AMDI < 0 & MPI > 0 ~ "M1-preactivated",
      AMDI > 0 & MPI > 0 ~ "M1-like",
      AMDI > 0 & MPI < 0 ~ "M2-like",
      TRUE ~ NA_character_ # for cases that do not match any condition
    ))

  myeloid_sub <- Seurat::SetIdent(myeloid_sub, value = myeloid_sub@meta.data$CellType)

  # Return to original object
  CellsMetaTrim <- subset(myeloid_sub@meta.data, select = c("MPI", "AMDI", "macrophage"))  #subset seurat object metadata; columns of interest

  main <- Seurat::AddMetaData(main, CellsMetaTrim)
  main$CellType_macspec <- factor(main$macrophage)

  for (msc_lvl in levels(main$CellType_macspec)) {
    cells_oi <- Seurat::Cells(subset(myeloid_sub, macrophage == msc_lvl))
    main <- Seurat::SetIdent(main, cells = cells_oi, value = paste(msc_lvl, "Macs"))
  }

  main$CellType_macspec <- Seurat::Idents(main)

  return(main)
}
