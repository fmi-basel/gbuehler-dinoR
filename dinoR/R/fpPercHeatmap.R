#' @title fpPercHeatmap
#'
#' @description Draws heatmaps of the percentages of all fragments in a ROI-sample combination in each footprint pattern.
#'
#' @details Draws heatmaps of the percentages of all fragments in a ROI-sample combination in each
#' footprint pattern supplied (for example: "tf", "open", "upNuc", "Nuc", "downNuc"). The rows of the heatmaps are split by ROI group.
#'
#' @param footprint_percentages A tibble where each column corresponds to a sample-footprint percentage and each row to a ROI,
#'  with the rows clustered by similarity.
#' @param breaks A list of vectors indicating numeric breaks used in (\code{ColorRamp2}) to define the heatmap color gradient,
#'  with one element per pattern (usually 5, or 3 if the nucleosome patterns have been combined).
#' @param plotcols A character vector of 5 colors to be used for the heatmaps of the 5 footprint
#' patterns ("tf", "open", "upNuc", "Nuc", "downNuc"), or 3 colors if the nucleosome patterns have been combined.
#'
#' @return Heatmaps of the percentages of all fragments in a ROI-sample combination in each footprint pattern.
#'
#' @examples
#' NomeData <- createExampleData()
#' NomeData <- footprintCalc(NomeData)
#' footprint_counts <- footprintQuant(NomeData)
#' fp <- footprintPerc(footprint_counts)
#' fpPercHeatmap(fp)
#'
#'
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom stringr str_extract
#'
#' @export
fpPercHeatmap <- function(footprint_percentages,breaks=rep(list(c(0,50,100)),5),
                          plotcols=c("#236467","#AA9B39","#822B56","#822B26","#822B99")){


    #patterns <- c("tf", "open", "upNuc", "Nuc", "downNuc")
    patterns <- unique(str_extract(colnames(footprint_percentages)[3:ncol(footprint_percentages)],"^[^_]+"))
    #initialize empty heatmap list
    ht_list <- NULL
    #draw heatmaps
    for (i in seq_along(patterns)){

        ht_list <- ht_list + Heatmap(matrix = as.matrix(footprint_percentages[,grep(paste0("^",patterns[i]),
                                                                              colnames(footprint_percentages))]),
                               border = TRUE,
                               col = colorRamp2(breaks = breaks[[i]],
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
