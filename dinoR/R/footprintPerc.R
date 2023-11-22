#' @title footprintPerc
#'
#' @description Calculates the percentage of all fragments in a ROI-sample combination corresponding to each footprint pattern.
#'
#' @details Calculates the percentage of all fragments in a ROI-sample combination in each footprint pattern. Then turns
#'  the table into wide format, where each column corresponds to a sample-footprint percentage and each row to a ROI
#'  and clusters the rows by similarity.
#'
#' @param footprint_counts A Summarized Experiment containing the sample names (colData), ROI names (rowData),
#' and number of fragments in each NOMe footprint pattern category as assays. For example
#' the output of the (\code{footprintQuant}) function.
#' @param minreads The minimum number of fragments to which a footprint could be assigned
#'  a ROI must have in all samples. All other ROIs are removed.
#' @param meanreads The minimum number of fragments to which a footprint could be assigned
#'  a ROI must have on average across all samples. All other ROIs are removed.
#' @param ROIgroup Column name of a metadata column in the (\code{rowData}) of the RSE,
#' describing a group each ROI belongs too, for example,
#' different transcription factor motifs at the center of the ROI.
#' @param combineNucCounts If TRUE, the upNuc, downNuc, and Nuc fragment counts will be combined into the Nuc category.
#'
#' @return A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI.
#'
#' @examples
#' NomeData <- createExampleData()
#' NomeData <- footprintCalc(NomeData)
#' footprint_counts <- footprintQuant(NomeData)
#' footprintPerc(footprint_counts)
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats dist hclust
#' @importFrom SummarizedExperiment assays assays<- rowData
#'
#' @export
footprintPerc <- function(footprint_counts, minreads=1,meanreads=1, ROIgroup="motif",combineNucCounts =FALSE){

  #keep only ROIs where at least one sample has more than 0 total counts
  footprint_counts <- footprint_counts[which(apply(assays(footprint_counts)[["all"]],1,min,na.rm=TRUE) >= minreads),]
  footprint_counts <- footprint_counts[which(apply(assays(footprint_counts)[["all"]],1,mean,na.rm=TRUE) >= meanreads),]

  if (combineNucCounts ==TRUE){
    # 3 patterns: calculate percentages  and re-arrange pattern quantification matrix from summarized experiment
    assays(footprint_counts)[["Nuc"]] <- assays(footprint_counts)[["Nuc"]] + assays(footprint_counts)[["upNuc"]] + assays(footprint_counts)[["downNuc"]]
    patterns <- c("tf","open","Nuc")
    assay_list <- list()
    for (p in seq_along(patterns)){
      #calculate percentages for 3 patterns
      assay_list[[p]] <- (assays(footprint_counts)[[patterns[p]]]/assays(footprint_counts)[["all"]])*100
      #remove the pattern from the column name
      colnames(assay_list[[p]]) <- paste(patterns[p], colnames(assay_list[[p]]),sep="_")
    }
  } else {
  # 5 patterns: calculate percentages  and re-arrange pattern quantification matrix from summarized experiment
  patterns <- c("tf","open","upNuc","Nuc","downNuc")
  assay_list <- list()
  for (p in seq_along(patterns)){
    #calculate percentages for all 5 patterns
    assay_list[[p]] <- (assays(footprint_counts)[[patterns[p]]]/assays(footprint_counts)[["all"]])*100
    #remove the pattern from the column name
    colnames(assay_list[[p]]) <- paste(patterns[p], colnames(assay_list[[p]]),sep="_")
   }
  }


  #combine into a data frame, together with ROI name and group
  patternQuantPercWide <- data.frame(ROI=rownames(rowData(footprint_counts)),
                                     ROIgroup=rowData(footprint_counts)[,ROIgroup],do.call("cbind",assay_list))

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
