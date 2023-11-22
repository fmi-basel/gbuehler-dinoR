#' @title metaPlots
#'
#' @description Plot the summarized GCH methylation protections across selected ROIs.
#'
#' @details  Summarizes the GCH methylation protections across selected ROIs.
#'
#' @param NomeData A Ranged Summarized Experiment (RSE) with an entry for each ROI. The (\code{rowData}) should contain information about each ROI,
#' including a ROIgroup.The (\code{assays}) should contain at least (\code{nFragsAnalyzed}) and (\code{reads}). (\code{nFragsAnalyzed}) describes the number of fragments
#' that were analyzed for each sample/ROI combination. (\code{reads}) contains a Gpos object for each sample/ROI combination,
#' with a position for each base in the ROI and two metadata columns (protection and methylation). protection is a sparse logical matrix where
#' TRUE stands for Cs protected from methylation, and methylation is a sparse logical matrix where TRUE stands for methylated Cs.
#' @param nr Integer used as a cutoff to filter sample ROI combinations that have less than
#' (\code{nr}) fragments analyzed (nFragsAnalyzed column).
#' @param nROI The number of ROIs that need to have a GpC methylation measurement at a given
#' position for this position to be included in the plot.
#' @param ROIgroup Column name of a metadata column in the rowData of the RSE,
#' describing a group each ROI belongs too, for example,
#' different transcription factor motifs at the center of the ROI.
#' @param span The (\code{span}) option to be used for the (\code{loess}) function (to draw a line through the datapoints).
#'
#' @return A tibble with the methylation protection profiles summarized across all ROIs in a certain group.
#'
#'
#' @examples
#' NomeData <- createExampleData()
#' metaPlots(NomeData)
#'
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom stats loess
#' @importFrom SummarizedExperiment colData rowData assays
#' @importFrom Matrix Matrix
#' @importFrom GenomicRanges mcols
#'
#' @export
metaPlots <- function(NomeData,nr=2,nROI=2,ROIgroup="motif",span=0.05){

  #get sample names and types and number of positions
  samples <- unique(colData(NomeData)$samples)
  types <- unique(rowData(NomeData)[,ROIgroup])
  npos <- length(assays(NomeData)[["reads"]][1,1][[1]])

  #average for each sample over all ROIs of a certain type
  perc.meth.list <- list()

  for (s in seq_along(samples)){
    #prepare empty matrix for average protection values...
    TypeAves <- matrix(nrow=length(types),ncol=npos)
    #...and number ROIs contributing to a certain position average
    TypeNs <- matrix(nrow=length(types),ncol=npos)

    # subset the RSE to the current sample
    NomeData1 <- NomeData[,samples[s]]

    for (t in seq_along(types)){
      #keep only the data for a given type of amplicon
      NomeData2 <-  NomeData1[rowData(NomeData1)[,ROIgroup] == types[t],]
      #filter out the amplicons that have < nr reads
      NomeData2 <-  NomeData2[as.vector(assays(NomeData2)[["nFragsAnalyzed"]] > nr),]

      #average protetction per position
      ampliAves <- vapply(assays(NomeData2)[["reads"]],.aveProt,rep(NaN,npos))

      #average over all the ROIs of a certain type t
      TypeAves[t,] <- rowMeans(ampliAves,na.rm=TRUE)*100
      TypeNs[t,] <- apply(X = ampliAves,1,function(x){length(x[!is.na(x)])})
    }
    #add names
    TypeAves <- data.frame(t(TypeAves))
    colnames(TypeAves) <- types

    TypeNs <- data.frame(t(TypeNs))
    colnames(TypeNs) <- types

    # add position headers
    ColNumbers <- -floor((nrow(TypeAves)/2)): floor((nrow(TypeAves)/2))
    if (nrow(TypeAves)%%2==0){
      TypeAves$position <- ColNumbers[-1]
      TypeNs$position <- ColNumbers[-1]
    } else {
      TypeAves$position <- ColNumbers
      TypeNs$position <- ColNumbers
    }

    #remove %protection values with less than n GpCs
    for (t in seq_along(types)){
      TypeAves[TypeNs[,t] < nROI,t] <- NaN
    }

    #make long
    perc.meth2 <- pivot_longer(TypeAves,cols=all_of(types),
                               names_to="type",values_to = "protection")
    perc.meth2 <- perc.meth2[is.na(perc.meth2$protection) ==FALSE,]

    #add smoothened values
    perc.meth <- data.frame()
    for (g in seq_along(unique(perc.meth2$type))){
      perc.meth3 <- perc.meth2[perc.meth2$type==unique(perc.meth2$type)[g],]
      perc.meth3$loess <- loess(perc.meth3$protection ~ perc.meth3$position,span=span)$fitted
      perc.meth <- rbind(perc.meth,perc.meth3)
    }
    #save the methylation profiles to a table
    perc.meth$sample <- samples[s]
    perc.meth.list[[s]] <- perc.meth
  }

  perc.meth <- tibble(do.call("rbind",perc.meth.list))
  return(perc.meth)
}

# function to combine the 2 sparse matrices (protection and methylation)...
# ...and calculate rowMeans (average protection per position)
.aveProt <- function(x){
  x_both <- mcols(x)[,"protection"] - mcols(x)[,"methylation"]
  x_both2 <- as.matrix(Matrix(x_both,sparse=FALSE))
  x_both2[x_both2==0] <- NA
  x_both2[x_both2==-1] <- 0
  rowMeans(x_both2,na.rm=TRUE)
}
