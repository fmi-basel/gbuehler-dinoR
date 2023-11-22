#' @title footprintQuant
#'
#' @description Quantify the footprint types.
#'
#' @details Count the number of fragments corresponding to a footprint type for each sample-ROI combination.
#'
#' @param NomeData A Ranged Summarized Experiment (RSE) with an entry for each ROI. The (\code{rowData}) should contain information about each ROI,
#' including a ROIgroup.The (\code{assays}) should contain at least (\code{nFragsAnalyzed}) and (\code{reads}). (\code{nFragsAnalyzed}) describes the number of fragments
#' that were analyzed for each sample/ROI combination. (\code{reads}) contains a Gpos object for each sample/ROI combination,
#' with a position for each base in the ROI and two metadata columns (protection and methylation). protection is a sparse logical matrix where
#' TRUE stands for Cs protected from methylation, and methylation is a sparse logical matrix where TRUE stands for methylated Cs.
#' In addition, there must be an assay called "footprints", which contains the assigned footprint ("tf","open","upNuc","Nuc","downNuc") for each fragment
#' (generated using the footprintCalc function).
#'
#' @return The Ranged Summarized Experiment with an assay added for each footprint type,
#'  containing the number of fragments that contain that footprint. An assay with the total number
#'  of pattern-able fragments ("all") is also added.
#' tf = transcription factor footprint
#' open = open chromatin footprint
#' upNuc = upstream nucleosome footprint
#' downNuc = downstream nucleosome footprint
#' Nuc = other nucleosome footprints
#'
#'
#' @examples
#' NomeData <- createExampleData()
#' NomeData <- footprintCalc(NomeData)
#' footprintQuant(NomeData)
#'
#' @importFrom SummarizedExperiment colData assays assay<-
#'
#' @export
footprintQuant <- function(NomeData){

  #extract samples
  samples <- unique(colData(NomeData)$samples)
  #define patterns
  patterns <- c("tf","open","upNuc","Nuc","downNuc")

  for (p in seq_along(patterns)){
  #calulcate patterns across ROIs for all samples, adding the patterns as assays
  assay(NomeData, patterns[p], withDimnames = FALSE) <- matrix(vapply(assays(NomeData)[["footprints"]],
                                                               function(x){length(which(x == patterns[p]))},
                                                               length(assays(NomeData)[["footprints"]])),ncol=length(samples))
  }

  #calulcate patternable reads across ROIs for all samples, adding the count as assay "all
  assay(NomeData, "all", withDimnames = FALSE) <- matrix(vapply(assays(NomeData)[["footprints"]],
                                                               function(x){length(which(!is.na(x)))},
                                                               length(assays(NomeData)[["footprints"]])),ncol=length(samples))

return(NomeData)

}
