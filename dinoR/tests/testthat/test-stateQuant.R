test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

library(tibble)
NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
  names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
  GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
  replace=TRUE),ncol=150,nrow=20)),20))
expect_no_error(stateQuant(NomeMatrix,nr=100))

NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
                     names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
                     GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=40*20,
                                                           replace=TRUE),ncol=40,nrow=20)),20))
expect_no_error(stateQuant(NomeMatrix,nr=2))
