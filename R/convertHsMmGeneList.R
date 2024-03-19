#' @title convertHsMmGeneList
#'
#' @description Converts between human and mouse gene names
#'
#' @param gene_list A string with 1 gene or a list with multiple genes
#' @param convert_to A string defining the direction of conversion (default is human)
#'
#' @return A list of genes
#'
#' @export
#'
convertHsMmGeneList <- function(gene_list, convert_to = "human"){
  human = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast", verbose = FALSE)
  mouse = biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast", verbose = FALSE)
  #human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  if (convert_to == "human") {
    genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene_list, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    return(humanx)
  } else {
    genesV2 = biomaRt::getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = gene_list, mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    mousex <- unique(genesV2[, 2])
    return(mousex)
  }
}
