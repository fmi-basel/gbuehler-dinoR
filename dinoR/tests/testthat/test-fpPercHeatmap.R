library(tibble)
NomeMatrix <- tibble(SampleName = c(rep("WT_1",5),
                                    rep("WT_2",5),rep("KO_1",5),rep("KO_2",5)),
                     names=rep(paste0("ROI",1:5),4),nFragsAnalyzed=rep(20,20),
                     GCH_DataMatrix=rep(c(list(matrix(rep(c(NA,0,1,1,NA,NA,NA,0,0,
                                                            1,NA,1,1,NA,1,NA,0,0,1,
                                                            NA,NA,0,0,1,NA,1,1,NA,0,0),
                                                          100),ncol=150,nrow=20,byrow=FALSE)),
                                          list(matrix(rep(c(NA,0,1,1,NA,1,NA,0,0,
                                                            1,NA,1,1,NA,1,NA,0,0,1,
                                                            NA,NA,0,0,1,NA,1,0,NA,0,0),
                                                          100),ncol=150,nrow=20,byrow=FALSE)),
                                          list(matrix(rep(c(NA,0,1,1,NA,NA,NA,0,0,
                                                            1,NA,1,1,1,1,NA,NA,1,1,
                                                            NA,NA,0,0,1,NA,1,1,1,0,0),
                                                          100),ncol=150,nrow=20,byrow=FALSE)),
                                          list(matrix(rep(c(1,0,NA,1,NA,NA,NA,0,0,
                                                            1,NA,0,0,NA,1,NA,0,0,1,
                                                            NA,NA,0,0,1,NA,1,1,NA,0,1),
                                                          100),ncol=150,nrow=20,byrow=FALSE)),
                                          list(matrix(rep(c(NA,0,1,1,0,0,NA,0,0,
                                                            1,NA,1,1,NA,1,NA,0,0,1,
                                                            NA,NA,0,0,1,NA,1,1,NA,0,1),
                                                          100),ncol=150,nrow=20,byrow=FALSE))),4))
footprint_counts <- footprintQuant(NomeMatrix)
fp <- footprintPerc(footprint_counts)

test_that("fpPercHeatmap works", {
  expect_no_error(fpPercHeatmap(fp))
})
