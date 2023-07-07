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
test_that("footprintPerc works", {
expect_no_error(footprintPerc(footprint_counts))
})

fp <- footprintPerc(footprint_counts)
fp <- tibble(cbind(fp[,1:2],round(fp[3:ncol(fp)],digits = 1)))
#make a tibble that resembles the output and compare
output <- tibble(ROI=c("ROI3","ROI1","ROI5","ROI2","ROI4"),
                 ROIgroup=rep("motif1",5),
                 tf_KO_1=rep(0,5),
                 tf_KO_2=rep(0,5),
                 tf_WT_1=rep(0,5),
                 tf_WT_2=rep(0,5),
                 open_KO_1=c(22.2,38.9,50.0,44.4,45.0),
                 open_KO_2=c(22.2,38.9,50.0,44.4,45.0),
                 open_WT_1=c(22.2,38.9,50.0,44.4,45.0),
                 open_WT_2=c(22.2,38.9,50.0,44.4,45.0),
                 upNuc_KO_1=rep(0,5),
                 upNuc_KO_2=rep(0,5),
                 upNuc_WT_1=rep(0,5),
                 upNuc_WT_2=rep(0,5),
                 Nuc_KO_1=c(72.2,50.0,44.4,50.0,50.0),
                 Nuc_KO_2=c(72.2,50.0,44.4,50.0,50.0),
                 Nuc_WT_1=c(72.2,50.0,44.4,50.0,50.00),
                 Nuc_WT_2=c(72.2,50.0,44.4,50.0,50.0),
                 downNuc_KO_1=c(5.6,11.1,5.6,5.6,5.0),
                 downNuc_KO_2=c(5.6,11.1,5.6,5.6,5.0),
                 downNuc_WT_1=c(5.6,11.1,5.6,5.6,5.0),
                 downNuc_WT_2=c(5.6,11.1,5.6,5.6,5.0)
                 )
test_that("footprintPerc returns the correct output", {
  expect_equal(fp,output)
})
