#' @title fpPercHeatmap
#'
#' @description Draws heatmaps of the percentages of all reads in a ROI-sample combination in each footprint pattern.
#'
#' @details Draws heatmaps of the percentages of all reads in a ROI-sample combination in each
#' footprint pattern ("tf", "open", "upNuc", "Nuc", "downNuc"). The rows of the heatmaps are split by ROI group.
#'
#' @param footprint_percentages A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI,
#'  with the rows clustered by similarity.
#' @param plotcols A character vector of 5 colors to be used for the heatmaps of the 5 footprint
#' patterns ("tf", "open", "upNuc", "Nuc", "downNuc").
#'
#' @return Heatmaps of the percentages of all reads in a ROI-sample combination in each footprint pattern.
#'
#' @examples
#' library(tibble)
#' counts <- tibble(sample = c(rep("WT",10),rep("KO",10)), ROI=rep(paste0("ROI",1:10),2),
#' tf=floor(runif(20,min=10,max=100)),open=floor(runif(20,min=10,max=100)),
#' upNuc=floor(runif(20,min=10,max=100)),Nuc=floor(runif(20,min=10,max=100)),
#' downNuc=floor(runif(20,min=10,max=100)))
#' randPerc <- footprintPerc(counts)
#' fpPercHeatmap(randPerc)
#'
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#'
#' @export
fpPercHeatmap <- function(footprint_percentages,
                          plotcols=c("#236467","#AA9B39","#822B56","#822B26","#822B99")){


    patterns <- c("tf", "open", "upNuc", "Nuc", "downNuc")
    #initialize empty heatmap list
    ht_list = NULL
    #draw heatmaps
    for (i in seq_along(patterns)){

        ht_list <- ht_list + Heatmap(matrix = as.matrix(footprint_percentages[,grep(paste0("^",patterns[i]),
                                                                              colnames(footprint_percentages))]),
                               border = TRUE,
                               col = colorRamp2(breaks = c(0,50,100),
                                                colors = c("white",plotcols[i],"black")),
                               cluster_rows = FALSE,
                               column_title = patterns[i],
                               name = patterns[i],
                               cluster_columns = FALSE,
                               show_column_names = TRUE,
                               show_row_names = FALSE,
                               na_col = "grey",
                               use_raster = FALSE,
                               heatmap_legend_param = list(title=patterns[i]),
                               split=footprint_percentages$ROIgroup
        )
    }
    return(ht_list)
}
