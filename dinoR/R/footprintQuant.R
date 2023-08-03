#' @title footprintQuant
#'
#' @description Count the number of reads in each NOMe state.
#'
#' @details Selects 3 windows (-50:-25, -8:8, 25:50) around the center of the provided region of interest (ROI) and calculates the average GpC
#' methylation protection for a given fragment across all GpCs in each window.
#' If it is above 0.5 the window is deemed protected, below 0.5, unprotected Depending on the protection pattern
#' in all windows, a read is put into one of 5 footprint categories: tf bound (0 - 1 - 0), open chromatin (0 - 0 - 0),
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
#' @param window_1 Integer vector with two elements representing start and end positions
#' of the first window relative to the ROI center.
#' @param window_2 Integer vector with two elements representing start and end positions
#' of the second window relative to the ROI center.
#' @param window_3 Integer vector with two elements representing start and end positions
#' of the third window relative to the ROI center.
#'
#' @return A SummarizedExperiment for the number of fragments in each footprint category,
#' with an assay for each pattern (see below), a column for each sample, and a row for each ROI.
#' tf = transcription factor footprint
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
#' footprintQuant(NomeMatrix)
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#'
#' @export
footprintQuant <- function(NomeMatrix,nr=2,window_1=c(-50,-25),window_2=c(-8,8),window_3=c(25,50)){
    NomeMatrix <- NomeMatrix[NomeMatrix$nFragsAnalyzed > nr,] #select only amplicon/sample pairs with >nr reads

    amplicons <- unique(NomeMatrix$names)
    if( length(amplicons) < 1){
        stop("Please provide ROI name(s)")
    }

    all.samples2 <- unique(NomeMatrix$SampleName)
    if( length(all.samples2) < 1){
        stop("Please provide sample name(s)")
    }

    #empty vectors for pattenrs
    quant1 <- character(0)
    quant2 <- character(0)
    quant3 <- numeric(0)
    quant4 <- numeric(0)
    quant5 <- numeric(0)
    quant6 <- numeric(0)
    quant7 <- numeric(0)
    quant8 <- numeric(0)

    for(am in seq_along(amplicons)){
        am1 <- NomeMatrix[NomeMatrix$names ==amplicons[am],]

        #extract available samples
        samples <- am1$SampleName

        #extract the list of 0/1 matrices
        peak1.list <- am1$GCH_DataMatrix

        #make sure there is a GCH matrix for every sample
        if( length(samples) != length(peak1.list)){
             stop(sprintf("Please provide a GCH protection matrix for every sample in ROI %s",amplicons[am]))
        }

        #add row.names indicating sample and make a corresponding vector (groups) indicating samples
        groups <- character(0)
        readIDs <- character(0)
        for (i in seq_along(peak1.list)){
            readIDs <- c(readIDs,row.names(peak1.list[[i]]))
            row.names(peak1.list[[i]]) <- paste(samples[i],seq(1,nrow(peak1.list[[i]])),sep="_")
            groups <- c(groups,rep(samples[i],nrow(peak1.list[[i]])))
        }

        #combine all matrices
        gch_protect <- do.call("rbind",peak1.list)

        #remove rows which have only NAs
        inx.row.NA <- which(rowMeans(gch_protect,na.rm=TRUE) != "NaN")
        gch_protect <- gch_protect[inx.row.NA,]
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

            #move it to a new object
            gch_protect2 <- gch_protect

            #check if we have at least 100 postions available in GCH prection matrix to cover all windows
            if(ncol(gch_protect2) < 100){
                warning(sprintf("GCH protection matrix for ROI %s will not cover all windows \n", amplicons[am]))
            }

            #label the columns with distance from center
            ColNumbers <- -floor((ncol(gch_protect2)/2)): floor((ncol(gch_protect2)/2))
            if (ncol(gch_protect2)%%2==0){
                colnames(gch_protect2) <- ColNumbers[-1]
            } else {
                colnames(gch_protect2) <- ColNumbers
            }

            #quantify 5 states: tf, open, nucleosome, uostream nucleosome, downstream nucleosome
            # quantify average methylation in 3 windows: default is [-50:-25], [-8:8], [25:50]

            #make a vector for annotating the selected windows
            windows <- ifelse(as.numeric(colnames(gch_protect2)) < window_1[1],"no window",
                        ifelse(as.numeric(colnames(gch_protect2)) >= window_1[1] &
                                 as.numeric(colnames(gch_protect2)) < window_1[2],"window1",
                               ifelse(as.numeric(colnames(gch_protect2)) >= window_1[2] &
                                        as.numeric(colnames(gch_protect2)) < window_2[1],"no window",
                                      ifelse(as.numeric(colnames(gch_protect2)) >= window_2[1] &
                                               as.numeric(colnames(gch_protect2)) <= window_2[2], "window2",
                                             ifelse(as.numeric(colnames(gch_protect2)) > window_2[2] &
                                                      as.numeric(colnames(gch_protect2)) <= window_3[1], "no window",
                                                    ifelse(as.numeric(colnames(gch_protect2)) > window_3[1] &
                                                             as.numeric(colnames(gch_protect2)) <= window_3[2], "window3",
                                                           "no window"))))))


           window1 <- round(rowMeans(gch_protect2[,windows=="window1"],na.rm=TRUE))
           window2 <- round(rowMeans(gch_protect2[,windows=="window2"],na.rm=TRUE))
           window3 <- round(rowMeans(gch_protect2[,windows=="window3"],na.rm=TRUE))

           #identify patterns
           #010=tf
           #000=accesible
           #all other patterns: nucleosome
           pattern <- ifelse(window1==0 & window2==1 & window3==0, "tf",
                        ifelse(window1==0 & window2==0 & window3==0, "open",
                               ifelse(window1==1 & window2==0 & window3==0, "Nuc",
                                      ifelse(window1==1 & window2==1 & window3==0, "upNuc",
                                             ifelse(window1==1 & window2==1 & window3==1, "Nuc",
                                                    ifelse(window1==0 & window2==1 & window3==1, "downNuc",
                                                           ifelse(window1==0 & window2==0 & window3==1, "Nuc",
                                                                  ifelse(window1==1 & window2==0 & window3==1, "Nuc",
                                                                         NA))))))))
           names(pattern) <- row.names(gch_protect2)

           #loop through samples
           for (i in seq_along(samples)){
               #select indexes of rows correspodning to current sample
               pattern3 <- pattern[which(groups == samples[i])]

               #quantify patterns
               quant1 <- c(quant1,amplicons[am])
               quant2 <- c(quant2,samples[i])
               quant3 <- c(quant3,length(which(pattern3 == "tf")))
               quant4 <- c(quant4,length(which(pattern3 == "open")))
               quant5 <- c(quant5,length(which(pattern3 == "Nuc")))
               quant6 <- c(quant6,length(which(is.na(pattern3))))
               quant7 <- c(quant7,length(which(pattern3 == "upNuc")))
               quant8 <- c(quant8,length(which(pattern3 == "downNuc")))
           }
        }
    }

    #save pattern quantifications as a summarized experiment...
    #... make a tibble with all data and make it wider to seperate the samples into columns
    pattern_quant <- data.frame(sample=quant2,ROI=quant1,tf=quant3,open=quant4,upNuc=quant7,Nuc=quant5,downNuc=quant8,noData=quant6)
    pattern_quant$all <- apply(pattern_quant[,3:7],1,sum)
    all_patterns <- colnames(pattern_quant)[3:9]
    pattern_quant_wide <- pivot_wider(pattern_quant,id_cols="ROI",names_from="sample",
                                      values_from=c("tf","open","upNuc","Nuc","downNuc","noData","all"),
                                      names_sort=TRUE,values_fill = 0)

    #....generate assay matrices for each pattern
    countsList <- list()
    for (p in seq_along(all_patterns)){
      countsList[[p]] <- as.matrix(pattern_quant_wide[,grep(paste0("^",all_patterns[p]),colnames(pattern_quant_wide))])
      colnames(countsList[[p]]) <- gsub(paste0(all_patterns[p],"_"),"",colnames(countsList[[p]]))
      rownames(countsList[[p]]) <- pattern_quant_wide$ROI
    }
    names(countsList) <- all_patterns

    #...generate summarized experiment
    sumExp <- SummarizedExperiment::SummarizedExperiment(assays=countsList,
                                                         rowData=pattern_quant_wide[,1], colData=data.frame(samples=colnames(countsList[[1]])))

    return(sumExp)
}
