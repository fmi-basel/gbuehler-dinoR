library(tibble)
library(SummarizedExperiment)

NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
  names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
  GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
  replace=TRUE),ncol=150,nrow=20)),20))

test_that("footprintQuant works", {
expect_no_error(footprintQuant(NomeMatrix,nr=2))
})
test_that("footprintQuant gives an error when nr > number of reads in all sample-ROI combinations", {
  expect_error(footprintQuant(NomeMatrix,nr=100))
})

NomeMatrix <- tibble(SampleName = c(rep("WT",10),rep("KO",10)),
                     names=rep(paste0("ROI",1:10),2),nFragsAnalyzed=rep(20,20),
                     GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=40*20,
                                                           replace=TRUE),ncol=40,nrow=20)),20))
test_that("footprintQuant warns that a region width is not covering all 3 windows", {
expect_warning(footprintQuant(NomeMatrix,nr=2),"will not cover all windows")
})

#check that the output is as expected
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

#construct expected output
#rowdata
rowdata <- data.frame(ROI=paste0("ROI",1:5))
rownames(rowdata) <- rowdata$ROI
#coldata
coldata <- data.frame(samples=c("KO_1","KO_2","WT_1","WT_2"))
rownames(coldata) <- coldata$samples
#assays
countsList <- list(
matrix(rep(0,20),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples)),
matrix(rep(c(7,8,4,9,9),4),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples)),
matrix(rep(0,20),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples)),
matrix(rep(c(9,9,13,10,8),4),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples)),
matrix(rep(c(2,1,1,1,1),4),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples)),
matrix(rep(0,20),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples)),
matrix(rep(c(18,18,18,20,18),4),ncol=4,nrow=5,dimnames = list(rowdata$ROI,coldata$samples))
)
names(countsList) <- c("tf", "open", "upNuc", "Nuc", "downNuc", "noData", "all")

output <- SummarizedExperiment(assays=countsList,
                     rowData=rowdata, colData=coldata)

test_that("footprintQuant returns the correct output", {
  expect_equal(footprint_counts,output)
})
