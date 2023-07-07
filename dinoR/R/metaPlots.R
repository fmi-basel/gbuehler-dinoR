#' @title metaPlots
#'
#' @description Plot the summarized GpC methylation protections across selected ROIs.
#'
#' @details  Summarizes the GpC methylation protections across selected ROIs.
#'
#' @param NomeMatrix A tibble where each row corresponds to a sample ROI combination,
#' which should contain the following columns:
#' nFragsAnalyzed = the number fragments (read pairs) that was available
#' names = the names of the ROIs that were analyzed
#' SampleName = the names of the samples that were analzyed
#' GCH_DataMatrix = for every combination of ROI and sample, a matrix where the columns
#' correspond to genomic positions from start to end of the ROI and the rows to fragments.
#' 1 = protected from GpC methylation. 0 = not protected from GpC methylation.
#' @param nr Integer used as a cutoff to filter sample ROI combinations that have less than
#' (\code{nr}) fragments analyzed (nFragsAnalyzed column).
#' @param nROI The number of ROIs that need to have a GpC methylation measuremnet at a given
#' position for this position to be included in the plot.
#' @param ROIgroup A vector of the same length as the number of rows in NomeMatrix,
#' describing a group each ROI belongs too, for example,
#' different transcription factor motifs at the center of the ROI.
#'
#' @return A tibble with the methylation protection profiles summarized across all ROIs in a certain group.
#'
#'
#' @examples
#' library(tibble)
#' NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
#' names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
#' GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
#' replace=TRUE),ncol=150,nrow=20)),20))
#' metaPlots(NomeMatrix)
#'
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect all_of
#' @importFrom stats loess
#'
#' @export
metaPlots <- function(NomeMatrix,nr=2,nROI=2,ROIgroup=rep("motif1",nrow(NomeMatrix))){
    NomeMatrix$type <- ROIgroup
    NomeMatrix <- NomeMatrix[NomeMatrix$nFragsAnalyzed > nr,] #select only amplicon/sample pairs with >nr reads

    amplicons <- unique(NomeMatrix$names)
    all.samples2 <- unique(NomeMatrix$SampleName)
    types <- unique(NomeMatrix$type)

    perc.meth.list <- list()

    #loop through samples
    for (s in seq_along(all.samples2)){
        Nome_matrix.shifted2 <- NomeMatrix[NomeMatrix$SampleName==all.samples2[s],]

        #loop through amplicon types
        TypeAves <- matrix(nrow=length(types),ncol=ncol(Nome_matrix.shifted2$GCH_DataMatrix[[1]]))
        TypeNs <- matrix(nrow=length(types),ncol=ncol(Nome_matrix.shifted2$GCH_DataMatrix[[1]]))

        for (t in seq_along(types)){
            NomeList <- Nome_matrix.shifted2$GCH_DataMatrix[Nome_matrix.shifted2$type==types[t]]
            #column means of all positions across reads per amplicon
            ampliAves <- t(vapply(NomeList,colMeans,rep(0,ncol(TypeAves)),na.rm=TRUE))
            #column means across amplicons for each position
            TypeAves[t,] <- colMeans(ampliAves,na.rm=TRUE)*100
            #number of amplicons with a 0/1 value (GpC) at each position
            TypeNs[t,] <- apply(X = ampliAves,2,function(x){length(x[!is.na(x)])})
       }

       #add names
       TypeAves <- data.frame(t(TypeAves))
       colnames(TypeAves) <- types

       TypeNs <- data.frame(t(TypeNs))
       colnames(TypeNs) <- types

       # add position headeers
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
       perc.meth2 <- pivot_longer(TypeAves,cols=all_of(types),names_to="type",values_to = "protection")
       perc.meth2 <- perc.meth2[is.na(perc.meth2$protection) ==FALSE,]

       perc.meth <- data.frame()
       for (g in seq_along(unique(perc.meth2$type))){
           perc.meth3 <- perc.meth2[perc.meth2$type==unique(perc.meth2$type)[g],]
          perc.meth3$loess <- loess(perc.meth3$protection ~ perc.meth3$position,span=0.05)$fitted
          perc.meth=rbind(perc.meth,perc.meth3)
       }

       #save the methylation profiles to a table
       perc.meth$sample <- all.samples2[s]
       perc.meth.list[[s]] <- perc.meth
    }

    perc.meth <- tibble(do.call("rbind",perc.meth.list))
    return(perc.meth)
}
