#' @title compareFootprints
#'
#' @description Compare each footprint pattern in WT and KO samples (percentages and diNOMeTest results).
#'
#' @details Plots the percentages of reads in each ROI in WT versus KO samples (mean of two replicates) in each
#' footprint pattern ("tf", "open", "upNuc", "Nuc", "downNuc"). The color indicates the ROI group and the shape the results of the diNOMeTest.
#'
#' @param footprint_percentages A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI,
#'  with the rows clustered by similarity.
#' @param res A tibble with the results of differential fragment count testing for each ROI-footprint combination.
#' @param WTsamples The control sample names (two replicates) as they appear in footprint_percentages.
#' @param KOsamples The treatment sample names (two replicates) as they appear in footprint_percentages.
#' @param plotcols A character vector of colors to be used for distinguishing the ROI groups (has to be the same length as there are ROI groups).
#'
#' @return A scatter plot for each footprint pattern comparing WT and KO percentages and significance test results.
#'
#' @examples
#' library(tibble)
#' NomeMatrix <- tibble(SampleName = c(rep("WT_1",5),
#' rep("WT_2",5),rep("KO_1",5),rep("KO_2",5)),
#' names=rep(paste0("ROI",1:5),4),nFragsAnalyzed=rep(20,20),
#' GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
#' replace=TRUE),ncol=150,nrow=20)),20))
#' footprint_counts <- footprintQuant(NomeMatrix)
#' res <- diNOMeTest(footprint_counts,WTsamples = c("WT_1","WT_2"),
#' KOsamples = c("KO_1","KO_2"))
#' footprint_percentages <- footprintPerc(footprint_counts)
#' compareFootprints(footprint_percentages,res,plotcols="black")
#'
#'
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom rlang .data
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_abline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 scale_shape_manual
#' @importFrom ggplot2 ggtitle
#' @importFrom cowplot plot_grid
#' @importFrom dplyr left_join
#'
#' @export
compareFootprints <- function(footprint_percentages, res, WTsamples = c("WT_1","WT_2"),
                              KOsamples = c("KO_1","KO_2"), plotcols){


  patterns <- c("tf", "open", "upNuc", "Nuc", "downNuc")

  plotlist <- list()
  for (i in seq_along(patterns)){
    #extract all columns with tf (or other pattern) starting
    patternQuantPercSel <- footprint_percentages[,grep(paste0("^",patterns[i]), colnames(footprint_percentages))]
    #remove tf_ to get sample names
    colnames(patternQuantPercSel) <- gsub(paste0(patterns[i],"_"),"",colnames(patternQuantPercSel))

    #average samples based on WT and KO specifications
    patternQuantPercSelAve <- data.frame(footprint_percentages[,1:2],
                                         WT=apply(patternQuantPercSel[,which(colnames(patternQuantPercSel) %in% WTsamples)],1,mean),
                                         KO=apply(patternQuantPercSel[,which(colnames(patternQuantPercSel) %in% KOsamples)],1,mean))

    # combine with tf_ contrast
    patternQuantPercSelAve <- left_join(patternQuantPercSelAve,res[res$contrasts==paste0(patterns[i],"_vs_all"),],by=c("ROI"="ROI"))

    #plot
    plotlist[[i]] <-  ggplot(patternQuantPercSelAve,aes(x=.data$WT,y=.data$KO,col=.data$ROIgroup,shape=.data$regulated)) +
      geom_point() + theme_classic() + xlim(c(0,100)) + ylim(c(0,100))+
      geom_abline(intercept=0,slope=1,linetype="dashed",col="grey",alpha=0.5) +
      scale_color_manual(values=plotcols) + scale_shape_manual(values=c("down"= 6,"no" = 0,"up" = 2)) + ggtitle(patterns[i])
  }

  return(plot_grid(plotlist = plotlist, ncol = 3, align = "vh"))
}

