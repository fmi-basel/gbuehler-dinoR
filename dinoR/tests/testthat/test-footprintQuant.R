library(SummarizedExperiment)

NomeData <- createExampleData()
NomeData <- footprintCalc(NomeData)

test_that("footprintQuant works", {
expect_no_error(footprintQuant(NomeData))
})

# check that the output is as expected
NomeData <- createExampleData(samples=c("WT_1"),group=c("WT"),nROI=10,randomMeth=FALSE)
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)

test_that("footprintQuant returns 9 assays", {
  expect_equal(length(assays(NomeData)),9)
})

test_that("footprintQuant returns 9 assays with correct names", {
  expect_equal(names(assays(NomeData)),c("nFragsAnalyzed","reads","footprints","tf","open","upNuc","Nuc","downNuc","all"))
})

TF_assay <- matrix(rep(5,10),ncol=1,nrow=10,dimnames=list(paste0("ROI",1:10),"WT_1"))

test_that("footprintQuant returns correct output in tf assay", {
  expect_equal(assays(NomeData)[["tf"]],TF_assay)
})
