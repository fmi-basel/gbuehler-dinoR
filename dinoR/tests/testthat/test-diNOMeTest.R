NomeData <- createExampleData()
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)

test_that("diNOMeTest works", {
  expect_no_error(diNOMeTest(NomeData))
})

# library(tibble)
# NomeMatrix <- tibble(SampleName = c(rep("WT_1",3),
#                                     rep("WT_2",3),rep("KO_1",3),rep("KO_2",3)),
#                      names=rep(paste0("ROI",1:3),4),nFragsAnalyzed=rep(20,12),
#                      GCH_DataMatrix=rep(c(list(matrix(rep(c(NA,0,1,1,NA,NA,NA,0,0,
#                                                             1,NA,1,1,NA,1,NA,0,0,1,
#                                                             NA,NA,0,0,1,NA,1,1,NA,0,0),
#                                                           100),ncol=150,nrow=20,byrow=FALSE)),
#                                           list(matrix(rep(c(NA,0,1,1,NA,1,NA,0,0,
#                                                             1,NA,1,1,NA,1,NA,0,0,1,
#                                                             NA,NA,0,0,1,NA,1,0,NA,0,0),
#                                                           100),ncol=150,nrow=20,byrow=FALSE)),
#                                           list(matrix(rep(c(NA,0,1,1,0,0,NA,0,0,
#                                                             1,NA,1,1,NA,1,NA,0,0,1,
#                                                             NA,NA,0,0,1,NA,1,1,NA,0,1),
#                                                           100),ncol=150,nrow=20,byrow=FALSE))),4))
# footprint_counts <- footprintQuant(NomeMatrix)
#
# test_that("diNOMeTest works", {
# expect_no_error(diNOMeTest(footprint_counts,WTsamples = c("WT_1","WT_2"),
# KOsamples = c("KO_1","KO_2")))
# })
#
# res <- diNOMeTest(footprint_counts,WTsamples = c("WT_1","WT_2"),
#            KOsamples = c("KO_1","KO_2"))
# res$logFC <- round(res$logFC)
# res$logCPM <- round(res$logCPM)
# res$F <- round(res$F)
#
# output <- tibble(
#   logFC= rep(0,15), logCPM = rep(16,15), F = rep(0,15),PValue=rep(1,15),
#   FDR=rep(1,15),contrasts=c(rep("open_vs_all",3),rep("tf_vs_all",3),rep("Nuc_vs_all",3),
#                                              rep("downNuc_vs_all",3),rep("upNuc_vs_all",3)),
# ROI=c("ROI2", "ROI3", "ROI1", "ROI1", "ROI2", "ROI3", "ROI3", "ROI2", "ROI1", "ROI1", "ROI2", "ROI3", "ROI1", "ROI2", "ROI3"),
# logadjPval=rep(0,15),regulated=rep("no",15)
# )
#
# test_that("diNOMeTest returns the correct output", {
#   expect_equal(res,output)
# })
