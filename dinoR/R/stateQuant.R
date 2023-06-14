#' @title stateQuant
#'
#' @description Count the number of reads in each NOMe state.
#'
#' @details Selects 3 windows (-50:-25, -8:8, 25:50) around the center of the provided region of interest (ROI) and calculates the average GpC
#' methylation protection for a given fragment across all GpCs in each window.
#' If it is above 0.5 the window is deemed protected, below 0.5, unprotected Depending on the protection pattern
#' in all windows, a read is put into one of 5 footprint categories: TF bound (0 - 1 - 0), open chromatin (0 - 0 - 0),
#' upstream nucleosome (1 - 1 - 0), other nucleosome (1 - 1 - 1, 1 - 0 - 0, 0 - 0 - 1, 1 - 0 - 1), and downstream nucleosome (0 - 1 - 1).
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
#'
#' @return A tibble with a row for each sample - ROI combination, and the number of fragments in each footprint category.
#' TF = transcription factor footprint
#' open = open chromatin footprint
#' upNuc = upstream nucleosome footprint
#' downNuc = downstream nucleosome footprint
#' Nuc = other nucleosome footprints
#' no Data = fragments where no clear footprint could be assigned due to lack of data (GpCs) in one of the windows
#'
#'
#' @examples
#' library(tibble)
#' NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
#' names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
#' GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
#' replace=TRUE),ncol=150,nrow=20)),20))
#' stateQuant(NomeMatrix)
#'
#' @importFrom tibble tibble
#'
#' @export
stateQuant <- function(NomeMatrix,nr=2){
  NomeMatrix <- NomeMatrix[NomeMatrix$nFragsAnalyzed > nr,] #select only amplicon/sample pairs with >nr reads

  amplicons <- unique(NomeMatrix$names)
  all.samples2 <- unique(NomeMatrix$SampleName)

  #empty vectors for pattenrs
  quant1 <- character(0)
  quant2 <- character(0)
  quant3 <- numeric(0)
  quant4 <- numeric(0)
  quant5 <- numeric(0)
  quant6 <- numeric(0)
  quant7 <- numeric(0)
  quant8 <- numeric(0)

  patterns.scores.df <- data.frame()

  for(am in seq_along(amplicons)){
    am1 <- NomeMatrix[NomeMatrix$names ==amplicons[am],]

    #extarct avalibale samples and order them by all.sampels
    samples <- am1$SampleName

    #extract the list of 0/1 matrices
    peak1.list <- am1$GCH_DataMatrix

    #add row.names indication sample and make a corresponding vector (groups) indication samples
    groups <- character(0)
    readIDs <- character(0)
    for (i in seq_along(peak1.list)){
      readIDs <- c(readIDs,row.names(peak1.list[[i]]))
      row.names(peak1.list[[i]]) <- paste(samples[i],1:nrow(peak1.list[[i]]),sep="_")
      groups <- c(groups,rep(samples[i],nrow(peak1.list[[i]])))
    }

    #combine all matrices
    gch_protect <- do.call("rbind",peak1.list)

    #remove rows which have only NAs
    inx.row.NA <- which(rowMeans(gch_protect,na.rm=TRUE) != "NaN")
    gch_protect <-gch_protect[inx.row.NA,]
    groups <- groups[inx.row.NA]
    readIDs <- readIDs[inx.row.NA]

    #extract the samples that were actually used, after removing rows with all NA
    samples <- unique(groups)

    #if there is less than 2 reads left in the gch_protect matrix, skip that amplicon
    if (is.null(dim(gch_protect))){
      next
    } else {

      #find GpCs
      GCinx <- which(colMeans(gch_protect,na.rm=TRUE) != "NaN")

      #for heatmap of all sites
      gch_protect2 <- gch_protect
      #colnames(gch_protect2) <- -311:311
      ColNumbers <- -floor((ncol(gch_protect2)/2)): floor((ncol(gch_protect2)/2))
      if (ncol(gch_protect2)%%2==0){
        colnames(gch_protect2) <- ColNumbers[-1]
      } else {
        colnames(gch_protect2) <- ColNumbers
      }


      #----------------------------------------------------------------------------------------------------------------
      #            quantify 3 states: TF, nucleosome, none as in https://www.sciencedirect.com/science/article/pii/S1097276520307930?via%3Dihub ####
      #----------------------------------------------------------------------------------------------------------------
      # quantify average methylation in 3 windows original: [-35:-25], [-7:7], [25:35]
      #current: [-50:-25], [-8:8], [25:50]

      #make a vector for annotating the selected windows
      windows <- ifelse(as.numeric(colnames(gch_protect2)) < -50,"no window",
                        ifelse(as.numeric(colnames(gch_protect2)) >= -50 & as.numeric(colnames(gch_protect2)) < -25,"window1",
                               ifelse(as.numeric(colnames(gch_protect2)) >= -25 & as.numeric(colnames(gch_protect2)) < -8,"no window",
                                      ifelse(as.numeric(colnames(gch_protect2)) >= -8 & as.numeric(colnames(gch_protect2)) <= 8, "window2",
                                             ifelse(as.numeric(colnames(gch_protect2)) > 8 & as.numeric(colnames(gch_protect2)) <= 25, "no window",
                                                    ifelse(as.numeric(colnames(gch_protect2)) > 25 & as.numeric(colnames(gch_protect2)) <= 50, "window3",
                                                           "no window"))))))


      window1 <- round(rowMeans(gch_protect2[,windows=="window1"],na.rm=TRUE))
      window2 <- round(rowMeans(gch_protect2[,windows=="window2"],na.rm=TRUE))
      window3 <- round(rowMeans(gch_protect2[,windows=="window3"],na.rm=TRUE))

      #identify patterns
      #010=TF
      #000=accesible
      #all other patterns: nucleosome
      pattern <- ifelse(window1==0 & window2==1 & window3==0, "TF",
                        ifelse(window1==0 & window2==0 & window3==0, "open",
                               ifelse(window1==1 & window2==0 & window3==0, "Nuc",
                                      ifelse(window1==1 & window2==1 & window3==0, "upNuc",
                                             ifelse(window1==1 & window2==1 & window3==1, "Nuc",
                                                    ifelse(window1==0 & window2==1 & window3==1, "downNuc",
                                                           ifelse(window1==0 & window2==0 & window3==1, "Nuc",
                                                                  ifelse(window1==1 & window2==0 & window3==1, "Nuc",
                                                                         NA))))))))
      #table(pattern)
      names(pattern) <- row.names(gch_protect2)


      #loop through samples
      #change this to starting with amplicon name for next time
      for (i in seq_along(samples)){
        #select indexes of rows correspodning to current sample
        inx <- which(groups == samples[i])
        gch_protect4 <- gch_protect2[inx,]
        pattern3 <- pattern[inx]
        #score3 <- score[inx]


        #quantify patterns
        quant1 <- c(quant1,amplicons[am])
        quant2 <- c(quant2,samples[i])
        quant3 <- c(quant3,length(which(pattern3 == "TF")))
        quant4 <- c(quant4,length(which(pattern3 == "open")))
        quant5 <- c(quant5,length(which(pattern3 == "Nuc")))
        quant6 <- c(quant6,length(which(is.na(pattern3))))
        quant7 <- c(quant7,length(which(pattern3 == "upNuc")))
        quant8 <- c(quant8,length(which(pattern3 == "downNuc")))
      }
    }
  }

  #save pattern quantifications
  pattern_quant <- tibble(sample=quant2,ROI=quant1,TF=quant3,open=quant4,upNuc=quant7,Nuc=quant5,downNuc=quant8,noData=quant6)
return(pattern_quant)

}
