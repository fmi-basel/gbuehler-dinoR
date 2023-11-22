#' @title compareFootprints
#'
#' @description Compare each footprint pattern in WT and KO samples (percentages and diNOMeTest results).
#'
#' @details Plots the percentages of reads in each ROI in WT versus KO samples (mean of two replicates) in each
#' footprint pattern. The color indicates the ROI group and the shape the results of the diNOMeTest.
#'
#' @param footprint_percentages A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI,
#'  with the rows clustered by similarity.
#' @param res A tibble with the results of differential fragment count testing for each ROI-footprint combination.
#' @param WTsamples The control sample names.
#' @param KOsamples The treatment sample names.
#' @param plotcols A character vector of colors to be used for distinguishing the ROI groups (has to be the same length as there are ROI groups).
#' @param facetROIgroup If TRUE, split the plots for each pattern by ROI group.
#'
#' @return A scatter plot for each footprint pattern comparing WT and KO percentages and significance test results.
#'
#' @examples
#' NomeData <- createExampleData()
#' NomeData <- footprintCalc(NomeData)
#' footprint_counts <- footprintQuant(NomeData)
#' res <- diNOMeTest(footprint_counts,WTsamples = c("WT_1","WT_2"),
#' KOsamples = c("KO_1","KO_2"))
#' footprint_percentages <- footprintPerc(footprint_counts)
#' compareFootprints(footprint_percentages,res,plotcols="black")
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_classic xlim ylim geom_abline scale_color_manual scale_shape_manual ggtitle facet_wrap vars
#' @importFrom rlang .data
#' @importFrom cowplot plot_grid
#' @importFrom dplyr left_join
#' @importFrom stringr str_extract
#'
#' @export
compareFootprints <- function(footprint_percentages, res, WTsamples = c("WT_1","WT_2"),
                              KOsamples = c("KO_1","KO_2"), plotcols, facetROIgroup = FALSE){


  #patterns <- c("tf", "open", "upNuc", "Nuc", "downNuc")
  patterns <- unique(str_extract(colnames(footprint_percentages)[3:ncol(footprint_percentages)],"^[^_]+"))

  plotlist <- list()
  for (i in seq_along(patterns)){
    #extract all columns with tf (or other pattern) starting
    patternQuantPercSel <- footprint_percentages[,grep(paste0("^",patterns[i]),
                                                       colnames(footprint_percentages))]
    #remove tf_ to get sample names
    colnames(patternQuantPercSel) <- gsub(paste0(patterns[i],"_"),"",colnames(patternQuantPercSel))

    #average samples based on WT and KO specifications
    patternQuantPercSelAve <- data.frame(footprint_percentages[,seq(1,2)],
                                         WT=apply(patternQuantPercSel[,which(colnames(patternQuantPercSel)
                                                                             %in% WTsamples)],1,mean),
                                         KO=apply(patternQuantPercSel[,which(colnames(patternQuantPercSel)
                                                                             %in% KOsamples)],1,mean))

    # combine with tf_ contrast
    patternQuantPercSelAve <- left_join(patternQuantPercSelAve,
                                        res[res$contrasts==paste0(patterns[i],"_vs_all"),],by=c("ROI"="ROI"))

    #plot
    plotlist[[i]] <-  ggplot(patternQuantPercSelAve,
                             aes(x=.data$WT,y=.data$KO,col=.data$ROIgroup,shape=.data$regulated)) +
      geom_point() + theme_classic() + xlim(c(0,100)) + ylim(c(0,100))+
      geom_abline(intercept=0,slope=1,linetype="dashed",col="grey",alpha=0.5) +
      scale_color_manual(values=plotcols) +
      scale_shape_manual(values=c("down"= 6,"no" = 0,"up" = 2)) + ggtitle(patterns[i])
   #split the plots by ROIgroup
    if(facetROIgroup == TRUE){
    plotlist[[i]] <- plotlist[[i]] + facet_wrap(vars(.data$ROIgroup))
    }
  }

  return(plot_grid(plotlist = plotlist, ncol = 3, align = "vh"))
}

