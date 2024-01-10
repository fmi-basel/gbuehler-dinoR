#' @title footprintCalc
#'
#' @description Assign a footprint type to each fragment based on GCH protection values in pre-defined windows.
#'
#' @details Selects 3 windows (default is -50:-25, -8:8, 25:50) around the center of the provided region of interest (ROI) and calculates the average GCH
#' methylation protection for a given fragment across all GCHs in each window.
#' If it is above 0.5 the window is deemed protected, below 0.5, unprotected. Depending on the protection pattern
#' in all windows, a read is put into one of 5 footprint categories: tf bound (0 - 1 - 0), open chromatin (0 - 0 - 0),
#' downstream positioned nucleosome (1 - 1 - 0), other nucleosome (1 - 1 - 1, 1 - 0 - 0, 0 - 0 - 1, 1 - 0 - 1), and upstream positioned nucleosome (0 - 1 - 1).
#'
#' @param NomeData A Ranged Summarized Experiment (RSE) with an entry for each ROI. The (\code{rowData}) should contain information about each ROI,
#' including a ROIgroup. The (\code{assays}) should contain at least (\code{nFragsAnalyzed}) and (\code{reads}). (\code{nFragsAnalyzed}) describes the number of fragments
#' that were analyzed for each sample/ROI combination. (\code{reads}) contains a Gpos object for each sample/ROI combination,
#' with a position for each base in the ROI and two metadata columns (protection and methylation). protection is a sparse logical matrix where
#' TRUE stands for Cs protected from methylation, and methylation is a sparse logical matrix where TRUE stands for methylated Cs.
#' @param window_1 Integer vector with two elements representing start and end positions
#' of the first window relative to the ROI center.
#' @param window_2 Integer vector with two elements representing start and end positions
#' of the second window relative to the ROI center.
#' @param window_3 Integer vector with two elements representing start and end positions
#' of the third window relative to the ROI center.
#'
#' @return The Ranged Summarized Experiment with an assay "footprints" added, which contains a footprint type assigned to each fragment.
#'
#' @examples
#' NomeData <- createExampleData()
#' footprintCalc(NomeData)
#'
#' @importFrom SummarizedExperiment assays assay<- colData
#' @importFrom Matrix Matrix
#' @importFrom BiocGenerics cbind
#' @importFrom GenomicRanges mcols
#'
#' @export
footprintCalc <- function(NomeData,window_1=c(-50,-25),window_2=c(-8,8),window_3=c(25,50)){

  #extract samples and length of regions
  samples <- unique(colData(NomeData)$samples)
  npos <- length(assays(NomeData)[["reads"]][1,1][[1]])

  #create a vector with distances from center
  DistToCenter <- -floor((npos/2)): floor((npos/2))
  if (npos%%2==0){
    DistToCenter <- DistToCenter[-1]
  }

  #make a vector for annotating the selected windows
  windows <- ifelse(DistToCenter < window_1[1],"no window",
                    ifelse(DistToCenter >= window_1[1] &
                             DistToCenter < window_1[2],"window1",
                           ifelse(DistToCenter >= window_1[2] &
                                    DistToCenter < window_2[1],"no window",
                                  ifelse(DistToCenter >= window_2[1] &
                                           DistToCenter <= window_2[2], "window2",
                                         ifelse(DistToCenter > window_2[2] &
                                                  DistToCenter <= window_3[1], "no window",
                                                ifelse(DistToCenter > window_3[1] &
                                                         DistToCenter <= window_3[2], "window3",
                                                       "no window"))))))



  #loop through the samples and calulcate patterns across ROIs for each sample, adding the patterns as assay named footprints
  SE_List <- list()
  for (s in seq_along(samples)){
      SE_List[[s]] <- NomeData[,samples[s]]
      assay(SE_List[[s]], "footprints", withDimnames = FALSE)  <- cbind(lapply(
        assays(SE_List[[s]])[["reads"]],.calcPatterns,windows=windows))
  }
  #combine the RSEs for all samples again
  NomeData <- do.call(BiocGenerics::cbind,SE_List)
  return(NomeData)
}



#function to calculate patterns
.calcPatterns <- function(sereads,windows){
  #combine sparse matrixes to protection matrix
  x <- sereads
  x_both <- mcols(x)[,"protection"] - mcols(x)[,"methylation"]
  x_both2 <- as.matrix(Matrix(x_both,sparse=FALSE))
  x_both2[x_both2==0] <- NA
  x_both2[x_both2==-1] <- 0

  if(dim(x_both2)[2] == 1) {
    #calculate average methylation within each window for each read
    Win1 <- round(mean(x_both2[windows=="window1"],na.rm=TRUE))
    Win2 <- round(mean(x_both2[windows=="window2"],na.rm=TRUE))
    Win3 <- round(mean(x_both2[windows=="window3"],na.rm=TRUE))
  } else {
    #calculate average methylation within each window for each read
    Win1 <- round(colMeans(x_both2[windows=="window1",],na.rm=TRUE))
    Win2 <- round(colMeans(x_both2[windows=="window2",],na.rm=TRUE))
    Win3 <- round(colMeans(x_both2[windows=="window3",],na.rm=TRUE))
  }
  #identify patterns
  #010=tf
  #000=accesible
  #all other patterns: nucleosome
  pattern <- ifelse(Win1==0 & Win2==1 & Win3==0, "tf",
                    ifelse(Win1==0 & Win2==0 & Win3==0, "open",
                           ifelse(Win1==1 & Win2==0 & Win3==0, "Nuc",
                                  ifelse(Win1==1 & Win2==1 & Win3==0, "downNuc",
                                         ifelse(Win1==1 & Win2==1 & Win3==1, "Nuc",
                                                ifelse(Win1==0 & Win2==1 & Win3==1, "upNuc",
                                                       ifelse(Win1==0 & Win2==0 & Win3==1, "Nuc",
                                                              ifelse(Win1==1 & Win2==0 & Win3==1, "Nuc",
                                                                     NA))))))))
  return(pattern)
}
