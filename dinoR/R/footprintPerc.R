#' @title footprintPerc
#'
#' @description Calculates the percentage of all reads in a ROI-sample combination in each footprint pattern.
#'
#' @details Calculates the percentage of all reads in a ROI-sample combination in each footprint pattern. Then turns
#'  the table into wide format, where each column corresponds to a sample-footprint percentage and each row to a ROI
#'  and clusters the rows by similarity.
#'
#' @param footprint_counts A Summarized Experiment containing the sample names (colData), ROI names (rowData),
#' and number of fragments in each NOMe footprint pattern category (assays). For example
#' the output of the footprintQuant function.
#' @param ROIgroup A vector of the same length as the number of rows in footprint_counts,
#' describing a group each ROI belongs too, for example,
#' different transcription factor motifs at the center of the ROI.
#'
#' @return A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI.
#'
#' @examples
#' library(tibble)
#' NomeMatrix <- tibble(SampleName = c(rep("WT_1",5),
#' rep("WT_2",5),rep("KO_1",5),rep("KO_2",5)),
#' names=rep(paste0("ROI",1:5),4),nFragsAnalyzed=rep(20,20),
#' GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
#' replace=TRUE),ncol=150,nrow=20)),20))
#' footprint_counts <- footprintQuant(NomeMatrix)
#' footprintPerc(footprint_counts)
#'
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom SummarizedExperiment assays
#' @importFrom SummarizedExperiment rowData
#'
#' @export
footprintPerc <- function(footprint_counts,
                          ROIgroup=rep("motif1",nrow(footprint_counts))){

  #keep only ROIs where all samples have more than 0 total counts
  footprint_counts <- footprint_counts[which(apply(assays(footprint_counts)[["all"]],1,min,na.rm=TRUE) > 0),]
  ROIgroup <- ROIgroup[which(apply(assays(footprint_counts)[["all"]],1,min,na.rm=TRUE) > 0)]

  #claculate percentages  and re-arrange pattern quantification matrix from summarized experiment
  patterns <- c("tf","open","upNuc","Nuc","downNuc")
  assay_list <- list()
  for (p in seq_along(patterns)){
    #calculate percentages for all 5 patterns
    assay_list[[p]] <- (SummarizedExperiment::assays(footprint_counts)[[patterns[p]]]/SummarizedExperiment::assays(footprint_counts)[["all"]])*100
    #remove the pattern from the column name
    colnames(assay_list[[p]]) <- paste(patterns[p], colnames(assay_list[[p]]),sep="_")
  }
  #combine into a data frame, together with ROI name and group
  patternQuantPercWide <- data.frame(ROI=SummarizedExperiment::rowData(footprint_counts)$ROI,ROIgroup=ROIgroup,do.call("cbind",assay_list))

  #change NaNs to NAs
  for (i in 3:ncol(patternQuantPercWide)){
    NaNinx <- which(patternQuantPercWide[,i]=="NaN")
    patternQuantPercWide[NaNinx,i] <- NA
  }

  #remove rows with NAs or NaNs in all samples
  patternQuantPercWide <- patternQuantPercWide[which(rowMeans(
    patternQuantPercWide[,3:ncol(patternQuantPercWide)],na.rm=TRUE) != "NaN"),]


  ####cluster the tables###
  dist.tf <- dist(patternQuantPercWide[,3:ncol(patternQuantPercWide)], method="euclidean")
  #in case there are NAs in the dissimilarity object, replace them with the mean distance
  # to prevent clustering from failing
  dist.tf_mean <- mean(dist.tf, na.rm=TRUE)
  dist.tf[is.na(dist.tf)] <- dist.tf_mean
  clust.tf <- hclust(dist.tf, method="ward.D2")
  #sort the whole data.frame  based on clustering
  patternQuantPercWideSorted <- tibble(patternQuantPercWide[clust.tf$order,])
  return(patternQuantPercWideSorted)
}
