NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
                     names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
                     GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=40*20,
                                                           replace=TRUE),ncol=40,nrow=20)),20))
test_that("metaPlots works", {
  expect_no_error(metaPlots(NomeMatrix))
})