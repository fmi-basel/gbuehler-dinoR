library(tibble)
NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
                     names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
                     GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=300*20,
                                                           replace=TRUE),ncol=300,nrow=20)),20))
test_that("metaPlots works", {
  expect_no_error(metaPlots(NomeMatrix,span=0.08))
})

NomeMatrix <- tibble(SampleName = c(rep("WT",10)),
                     names=rep(paste0("ROI",1:10),1),nFragsAnalyzed=rep(20,10),
                     GCH_DataMatrix=rep(list(matrix(rep(c(NA,0,1,1,NA,NA,NA,0,0,
                                                          1,NA,1,1,NA,1,NA,0,0,1,NA),
                                                        40),ncol=40,nrow=20,byrow=FALSE)),10))
plotData <- metaPlots(NomeMatrix,span=0.2)
plotData$protection <- round(plotData$protection)
plotData$loess <- round(plotData$loess)
output <- tibble(position = -19:20,type=rep("motif1",40),
                 protection=rep(58,40),loess=rep(58,40),
                 sample=rep("WT",40))

test_that("metaPlots returns the correct output", {
  expect_equal(plotData,output)
})

