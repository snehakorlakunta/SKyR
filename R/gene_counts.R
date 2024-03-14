#' @title gene_counts
#'
#' @description Creates a data.frame of cell counts for a list of genes.
#'
#' @param seurat_object A Seurat object
#' @param gene_list A string with 1 gene or a list with multiple genes
#' @param group Metadata column of interest (default is active.ident)
#'
#' @return A data.frame with gene names, metadata, average expression, percent
#'     expression, and cell counts.
#'
#' @export
#'
gene_counts <- function(seurat_object, gene_list, group = "ident") {
  gene_list2 <- intersect(gene_list, rownames(seurat_object))
  if (length(gene_list) != length(gene_list2)) {
    warning("The following genes were not found in the RNA assay: ", setdiff(gene_list, gene_list2))
  }
  res <- NULL

  for (g in gene_list2) {
    avgexpr <- Seurat::AverageExpression(seurat_object, assay = "RNA", features = g, group.by = group)
    avgexpr <- avgexpr$RNA
    clean_colnames <- colnames(avgexpr)

    pctexpr <- scCustomize::Percent_Expressing(seurat_object, assay = "RNA", features = g, group_by = group)
    colnames(pctexpr) <- clean_colnames

    cellcts <- Seurat::FetchData(seurat_object, vars = c(group, g))
    cellcts <- cellcts[cellcts[[g]] > 0,] #subset(cellcts, Pdgfra > 0)
    cellcts <- table(cellcts[[group]])
    cellcts <- t(data.matrix(cellcts))
    if (sum(Seurat::GetAssayData(seurat_object, assay = "RNA")[g,]>0) != sum(cellcts)) {
      warning("Error in cell counts for ", g)
    }

    rownames(pctexpr) <- "Percent Expr"
    rownames(avgexpr) <- "Average Expr"
    rownames(cellcts) <- "N"

    pctexpr <- data.frame(pctexpr)
    avgexpr <- data.frame(avgexpr)
    cellcts <- data.frame(cellcts)

    tbl <- rbind(avgexpr, pctexpr)
    tbl <- rbind(tbl, cellcts)
    tbl <- data.frame(tbl)
    colnames(tbl) <- clean_colnames

    tbl["Feature",] <- g

    tbl <- t(tbl)
    tbl <- data.frame(tbl)

    tbl <- rownames_to_column(tbl, group)

    tbl <- tbl[, c(5, 1, 2, 3, 4)]

    if (is.null(res)) {
      res <- tbl
    } else {
      res <- rbind(res, tbl)
    }
  }
  return (res)
}
