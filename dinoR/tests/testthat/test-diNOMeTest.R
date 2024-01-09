library(tibble)

# check that it works with randomly generated data
NomeData <- createExampleData()
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)

test_that("diNOMeTest works", {
  expect_no_error(diNOMeTest(NomeData))
})

# check that the output is as expected withd default example data
NomeData <- createExampleData(samples=c("WT_1","WT_2","KO_1","KO_2"),group=c("WT","WT","KO","KO"),nROI=10,randomMeth=FALSE)
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)


res <- diNOMeTest(NomeData,WTsamples = c("WT_1","WT_2"),
            KOsamples = c("KO_1","KO_2"))
 res$logFC <- round(res$logFC)
 res$logCPM <- round(res$logCPM)
 res$F <- round(res$F)

 output <- tibble(
   logFC= rep(0,50), logCPM = rep(14,50), F = rep(0,50),PValue=rep(1,50),
   FDR=rep(1,50),contrasts=c(rep("open_vs_all",10),rep("tf_vs_all",10),rep("Nuc_vs_all",10),
                                              rep("downNuc_vs_all",10),rep("upNuc_vs_all",10)),
 ROI=rep(paste0("ROI",1:10),5),
 motif=rep("motif1",50),
 logadjPval=rep(0,50),regulated=rep("no",50)
 )

test_that("diNOMeTest returns the correct output", {
   expect_equal(res,output)
})
