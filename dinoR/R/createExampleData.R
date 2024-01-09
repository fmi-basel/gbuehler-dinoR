#' @title createExampleData
#'
#' @description Creates an RSE object with mock NOMe-seq data.
#'
#' @details Creates an RSE object with mock NOMe-seq data.
#'
#' @param samples The sample names.
#' @param group The sample group names.
#' @param nROI The number of ROIs that should be constructed.
#' @param randomMeth Logical indicating whether the methylation/protection values should be randomly generated.
#'
#' @return RSE object with mock data.
#'
#' @examples
#' createExampleData()
#'
#' @importFrom GenomicRanges GRanges GPos seqnames start end strand
#' @importFrom SummarizedExperiment SummarizedExperiment colData rowRanges assays assay<-
#' @importFrom Matrix Matrix
#'
#' @export
createExampleData <- function(samples=c("WT_1","WT_2","KO_1","KO_2"),group=c("WT","WT","KO","KO"),nROI=20,randomMeth=TRUE){

#construct a Ranged Summarized Experiment containing the NOMeseq data
#sample annotations (colData)
annots <- data.frame(samples=samples,group=group)
rownames(annots) <- annots$samples

#ROI annotations (rowData)
ROIs_gr <- rep(GRanges("chr1:100-400:+"),nROI)
ROIs_gr$motif <- "motif1"
names(ROIs_gr) <- paste0("ROI",seq_len(nROI))

#assay of number of frgamnets analyzed
assay_types <- c("nFragsAnalyzed")
assay_list <- list()
for (i in seq_along(assay_types)){
  NomeMatrix_wide <- matrix(ncol=length(samples),nrow=nROI,20)
  colnames(NomeMatrix_wide) <- samples
  rownames(NomeMatrix_wide) <- paste0("ROI",seq_len(nROI))
  assay_list[[i]] <- NomeMatrix_wide
}
names(assay_list) <- assay_types

#combine into RSE
NomeData <- SummarizedExperiment(colData = annots,
                                 rowRanges = ROIs_gr,
                                 assays = assay_list)

#add methylation info
SE_List <- list()
for(s in seq_along(annots$samples)){
  gr_list1 <- list()
  for (r in seq_along(ROIs_gr)){
    gpos1 <- GPos(seqnames=seqnames(ROIs_gr)[r], pos=start(ROIs_gr)[r]:end(ROIs_gr)[r],
                  strand=strand(ROIs_gr)[r],
                  seqinfo=NULL, seqlengths=NULL, stitch=NA)

    if(randomMeth == TRUE){ #randomly generate the TRUE FALSE values for the sparse matrices
           gpos1$protection <- Matrix(matrix(sample(c(TRUE,FALSE),size=301*20,
                                             replace=TRUE),nrow=301,ncol=20),sparse=TRUE)
           gpos1$methylation <- Matrix(matrix(sample(c(TRUE,FALSE),size=301*20,
                                              replace=TRUE),nrow=301,ncol=20),sparse=TRUE)
    } else { # keep the TRUE FALSE values constant (for testing)
           meth <- rep(c(TRUE,FALSE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE,FALSE,FALSE),301)
           prot <- rep(c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE),301)
           gpos1$protection <- Matrix(matrix(prot,nrow=301,ncol=20),sparse=TRUE)
           gpos1$methylation <- Matrix(matrix(meth,nrow=301,ncol=20),sparse=TRUE)
    }
    gr_list1[[r]] <- gpos1
  }
  names(gr_list1) <- names(ROIs_gr)
  SE_List[[s]] <- NomeData[,annots$samples[s]]
  assay(SE_List[[s]], "reads", withDimnames = FALSE)  <- cbind(gr_list1)
}
NomeData <- do.call("cbind",SE_List)
return(NomeData)
}
