#' @title footprintPerc
#'
#' @description Calculates the percentage of all reads in a ROI-sample combination in each footprint pattern.
#'
#' @details Calculates the percentage of all reads in a ROI-sample combination in each footprint pattern. Then turns
#'  the table into wide format, where each column corresponds to a sample-footprint percentage and each row to a ROI
#'  and clusters the rows by similarity.
#'
#'
#' @param footprint_counts A tibble containing the sample names, ROI names,
#' and number of fragments in each NOMe footprint pattern category. For example
#' the output of the stateQuant function. Rows without any footprint counts should be removed.
#' @param ROIgroup A vector of the same length as the number of rows in footprint_counts,
#' describing a group each ROI belongs too, for example,
#' different transcription factor motifs at the center of the ROI.
#'
#' @return A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI.
#'
#' @examples
#' library(tibble)
#' counts <- tibble(sample = c(rep("WT",10),rep("KO",10)), ROI=rep(paste0("ROI",1:10),2),
#' TF=floor(runif(20,min=10,max=100)),open=floor(runif(20,min=10,max=100)),
#' upNuc=floor(runif(20,min=10,max=100)),Nuc=floor(runif(20,min=10,max=100)),
#' downNuc=floor(runif(20,min=10,max=100)))
#' footprintPerc(counts)
#'
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom stats dist
#' @importFrom stats hclust
#'
#' @export
footprintPerc <- function(footprint_counts,
                          ROIgroup=rep("motif1",nrow(footprint_counts))){

  patternQuant <- footprint_counts[,1:7]
  #sum up all patterned reads
  patternQuant$total <- rowSums(patternQuant[3:7])

  #tell users to remove sample-ROI combinations with 0 total reads (all reads ended up in noData)
  if (length(which(patternQuant$total == 0)) > 0){
    stop("Please remove sample-ROI combination with zero counts from your input!")
  }

  #calculate %TF and % open and % nucleosome
  patternQuantPerc <- cbind(patternQuant[,1:2],
                            ROIgroup,(patternQuant[,3:7]/patternQuant[,c(8,8,8,8,8)])*100)

  #put each sample in a column
  patternQuantPercWide <- pivot_wider(patternQuantPerc,id_cols = c("ROI","ROIgroup"),
                                      names_from = "sample",values_from =
                                        c("TF", "open", "upNuc", "Nuc", "downNuc"),names_sort=TRUE)

  #change NaNs to NAs
  for (i in 3:ncol(patternQuantPercWide)){
    NaNinx <- which(patternQuantPercWide[,i]=="NaN")
    patternQuantPercWide[NaNinx,i] <- NA
  }

  #remove rows with NAs or NaNs in all samples
  patternQuantPercWide <- patternQuantPercWide[which(rowMeans(
    patternQuantPercWide[,3:ncol(patternQuantPercWide)],na.rm=TRUE) != "NaN"),]


  ####cluster the tables###
  dist.TF <- dist(patternQuantPercWide[,3:ncol(patternQuantPercWide)], method="euclidean")
  #in case there are NAs in the dissimilarity object, replace them with the mean distance
  # to prevent clustering from failing
  dist.TF_mean <- mean(dist.TF, na.rm=T)
  dist.TF[is.na(dist.TF)] <- dist.TF_mean
  clust.TF <- hclust(dist.TF, method="ward.D2")
  #sort the whole data.frame  based on clustering
  patternQuantPercWideSorted <- tibble(patternQuantPercWide[clust.TF$order,])
  return(patternQuantPercWideSorted)



}
