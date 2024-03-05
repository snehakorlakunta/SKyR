#' @title save_filt_seu
#'
#' @description Creates a filtered Seurat object
#'
#' @param seurat_obj Unfiltered Seurat object
#' @param counts_file Unprocessed barcode-feature-matrix file
#' @param save_file_name String for filename of saved RDS file
#' @param a_MIN_GENES_PER_CELL Lower threshold for genes per cell
#' @param a_MAX_GENES_PER_CELL Upper threshold for genes per cell
#' @param a_MIN_COUNTS_PER_CELL Lower threshold for counts per cell
#' @param a_MAX_COUNTS_PER_CELL Upper threshold for counts per cell
#' @param a_MAX_PCT_MITO Upper threshold for mitochondrial gene expression
#' @param a_MAX_PCT_RIBO Upper threshold for ribosomal gene expression
#'
#' @return A Seurat object with filtered values
#'
#' @export
#'
save_filt_seu <- function(seurat_obj, counts_file, save_file_name,
                          a_MIN_GENES_PER_CELL, a_MAX_GENES_PER_CELL,
                          a_MIN_COUNTS_PER_CELL, a_MAX_COUNTS_PER_CELL,
                          a_MAX_PCT_MITO, a_MAX_PCT_RIBO#, a_MAX_PCT_HGB, a_MAX_PCT_PLT
) {
  if(length(grep("^mt-", rownames(counts_file), value=T)) == 0) {
    if(length(grep("^MT-", rownames(counts_file), value=T)) == 0) {
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^Mt-", col.name = "percent_mito")           # mitochondria
    } else {
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^MT-", col.name = "percent_mito")           # mitochondria
    }
  } else {
    seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^mt-", col.name = "percent_mito")           # mitochondria
  }

  if(length(grep("^rp[sl]", rownames(counts_file), value=T)) == 0) {
    if(length(grep("^RP[SL]", rownames(counts_file), value=T)) == 0) {
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^Rp[sl]", col.name = "percent_ribo")        # ribosomes
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^Hb[^(p)]", col.name = "percent_hb")        # hemoglobin
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^Pecam1$|^Pf4$", col.name = "percent_plat")    # platelets
    } else {
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^RP[SL]", col.name = "percent_ribo")        # ribosomes
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^HB[^(P)]", col.name = "percent_hb")        # hemoglobin
      seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^PECAM1$|^PF4$", col.name = "percent_plat")    # platelets
    }
  } else {
    seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^rp[sl]", col.name = "percent_ribo")        # ribosomes
    seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^hb[^(p)]", col.name = "percent_hb")        # hemoglobin
    seurat_obj <- Seurat::PercentageFeatureSet(seurat_obj, "^pecam1$|^pf4$", col.name = "percent_plat")    # platelets
  }

  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > a_MIN_GENES_PER_CELL
                       & nFeature_RNA < a_MAX_GENES_PER_CELL
                       & nCount_RNA > a_MIN_COUNTS_PER_CELL
                       & nCount_RNA < a_MAX_COUNTS_PER_CELL
                       & percent_mito < a_MAX_PCT_MITO
                       & percent_ribo < a_MAX_PCT_RIBO
                       #& percent_hb < a_MAX_PCT_HGB
                       #& percent_plat < a_MAX_PCT_PLT
  )

  saveRDS(seurat_obj, paste0(save_file_name, ".filt.rds"))

  return (seurat_obj)
}
