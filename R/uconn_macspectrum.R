#' @title uconn_macspectrum
#'
#' @description Calculating MPI and AMDI of macrophage samples using MacSpectrum's algorithm.
#'
#' @param mac_mtx A data.frame containing the gene expression levels of your dataset.
#'     In the data.frame, rows are genes and columns are samples (e.g. individual single cells,
#'     regular RNA-seq samples, qPCR samples, etc. see example below). The first column contains
#'     Ensembl IDs of the genes, and other columns are gene expression levels of each cell or sample,
#'     in UMI, FPKM/RPKM, TPM format or other relative expression formats (eg. qPCR results).
#'     Currently, mouse or human Ensembl IDs are supported. Human genes in your sample dataset will be
#'     automatically mapped to their murine homologs using MGI homology. The colname of the first column
#'     should be "geneid", and colnames of the rest colnames should be sample names.
#' @param feature a character vector which contains any features of the samples (treated/control,
#'     disease/healthy, etc.). The features should be of the same number and order as the samples in
#'     the mac_mtx (e.g., cell1, cell2, cell3, etc.; or sample1, sample2, sample3, etc.).
#' @param select_hu_mo a string indicating if the samples are of human ("hum") or mouse ("mou"). default "mou"
#'
#' @references https://macspectrum.uconn.edu/
#'
#' @return a data.frame containing sample names, features, MPI, and AMDI
#'
#' @keywords internal
#' @export
uconn_macspectrum <- function(mac_mtx, feature, select_hu_mo="mou") { #select_hu_mo is either hum or mou

  #data("m.M1_mean")      #load(paste0(macspec_path, "data/m.M1_mean.rda"))
  #data("m.M2_mean")      #load(paste0(macspec_path, "data/m.M2_mean.rda"))
  #data("m.M0_mean")      #load(paste0(macspec_path, "data/m.M0_mean.rda"))
  #data("m.hum_mou_map")  #load(paste0(macspec_path, "data/m.hum_mou_map.rda"))

  M1_mean<-m.M1_mean
  M2_mean<-m.M2_mean
  M0_mean<-m.M0_mean
  hum_mou_map<-m.hum_mou_map

  rownames(M1_mean)<-M1_mean$GeneID
  rownames(M2_mean)<-M2_mean$GeneID
  rownames(M0_mean)<-M0_mean$GeneID

  inFile <- mac_mtx
  if (is.null(inFile)) {return(NULL)}
  #    inFile <-read.csv(inFile$datapath, header = header)
  #
  if (select_hu_mo=="hum"){
    inFile_gene_id<-1:nrow(inFile)
    for (i in 1:nrow(inFile)){
      inFile_gene_id[i]<-as.character((hum_mou_map[,2])[inFile[i,1] == hum_mou_map[,1]])[1]
    }
    inFile_gene_id[is.na(inFile_gene_id)]<-"no_match"
    inFile[,1]<-inFile_gene_id
  }

  inFile<-inFile[!duplicated(inFile[,1]),]#remove duplicated genes
  rownames(inFile)<-inFile[,1]

  inFile_feature<-feature
  if (is.null(inFile_feature)) {return(NULL)}
  #    inFile_feature <-read.csv(inFile_feature$datapath, header = header,stringsAsFactors=F)



  inFile<-inFile[,2:ncol(inFile)]#drop first column of gene ID
  inFile<-inFile-rowMeans(inFile)




  MPI_genes<-intersect(M1_mean$GeneID,rownames(inFile))
  M1_mean<-M1_mean[MPI_genes,]
  M2_mean<-M2_mean[MPI_genes,]
  inFile_bak<-inFile
  inFile<-inFile[MPI_genes,]

  AMDI_genes<-intersect(M0_mean$GeneID,rownames(inFile_bak))
  M0_mean<-M0_mean[AMDI_genes,]
  inFile_m0<-inFile_bak[AMDI_genes,]

  #sigma of mac cells:
  inFile_sigma<-1:ncol(inFile)
  total_gene_number<-nrow(inFile)
  for(i in 1:ncol(inFile)) {
    options(digits=9)
    inFile_sigma[i]<-(sum(inFile[,i]^2)/total_gene_number)^0.5
  }

  #sigma of mac cells for m0:
  inFile_sigma_m0<-1:ncol(inFile_m0)
  total_gene_number<-nrow(inFile_m0)
  for(i in 1:ncol(inFile_m0)) {
    options(digits=9)
    inFile_sigma_m0[i]<-(sum(inFile_m0[,i]^2)/total_gene_number)^0.5
  }

  #sigma of M0 mean:
  total_gene_number<-nrow(M0_mean)
  M0_sigma<-(sum(M0_mean$value^2)/total_gene_number)^0.5

  #sigma of M1 mean:
  total_gene_number<-nrow(M1_mean)
  M1_sigma<-(sum(M1_mean$value^2)/total_gene_number)^0.5

  #sigma of M2 mean:
  total_gene_number<-nrow(M2_mean)
  M2_sigma<-(sum(M2_mean$value^2)/total_gene_number)^0.5

  #correlation of Mac  to M0 mean:
  total_gene_number<-nrow(inFile_m0)
  inFile_Pearson_per_cell_m0<-1:ncol(inFile_m0)
  for (j in 1:ncol(inFile_m0)){
    inFile_Pearson_per_cell_m0[j]<-sum((inFile_m0[,j]/inFile_sigma_m0[j])*(M0_mean$value/M0_sigma))/total_gene_number
  }

  #correlation of Mac  to M1 mean:
  total_gene_number<-nrow(M2_mean)
  inFile_Pearson_per_cell_m1<-1:ncol(inFile)
  for (j in 1:ncol(inFile)){
    inFile_Pearson_per_cell_m1[j]<-sum((inFile[,j]/inFile_sigma[j])*(M1_mean$value/M1_sigma))/total_gene_number
  }

  #correlation of ATM  to M2 mean:
  total_gene_number<-nrow(M2_mean)
  inFile_Pearson_per_cell_m2<-1:ncol(inFile)
  for (j in 1:ncol(inFile)){
    inFile_Pearson_per_cell_m2[j]<-sum((inFile[,j]/inFile_sigma[j])*(M2_mean$value/M2_sigma))/total_gene_number
  }

  a<-0.991414467
  b<-1
  c<- -0.0185412856
  x0<-inFile_Pearson_per_cell_m1
  y0<-inFile_Pearson_per_cell_m2
  d_sqr<-(a*x0+b*y0+c)^2/(a^2+b^2)
  x_start<--1
  y_start<-(-a)*x_start+(-c)
  x_end<-1
  y_end<-(-a)*x_end+(-c)

  l<-((x0-x_start)^2+(y0-y_start)^2-d_sqr)^0.5
  l_max<-((x_end-x_start)^2+(y_end-y_start)^2-d_sqr)^0.5
  MPI<-(l-0)/(l_max-0)*100-50

  AMDI<- -inFile_Pearson_per_cell_m0*50

  mac_output<-data.frame(colnames(inFile),inFile_feature,MPI,AMDI,row.names=colnames(inFile),stringsAsFactors=F)
  colnames(mac_output)<-c("Samples","Feature","MPI","AMDI")
  return(mac_output)
}
