#' @title pct_gene_cat
#'
#' @description Creates a scatter plot of genes matching a particular
#'    pattern and a set threshold for max percent
#'
#' @param gene_cat String category for selected gene subset
#' @param pattern String for grep-ing genes
#' @param max_pct Threshold for max-cutoff for genes
#' @param print_genes Choose whether to print applicable genes, default is FALSE
#'
#' @return A sorted scatterplot with a set max threshold
#'
#' @export
pct_gene_cat <- function(gene_cat, pattern, max_pct, print_genes=FALSE) {
  cat_genes <- grep(pattern, rownames(counts_file), ignore.case=T, value=T)
  if (print_genes) print(cat_genes)

  # compute pct genes in category
  cat_gene_read_counts = Matrix::colSums(counts_file[cat_genes,])
  pct_cat = cat_gene_read_counts / counts_per_cell * 100
  plot(sort(pct_cat), main= paste0(dataset_name, ": ", gene_cat),
       xlab = paste0("cells sorted by percentage ", gene_cat, " counts"),
       ylab = paste0("percentage ", gene_cat, " counts"))

  abline(h=max_pct, col='red')
}
