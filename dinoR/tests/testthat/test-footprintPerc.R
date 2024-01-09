library(tibble)

#test that footprintPerc works with random example data
NomeData <- createExampleData()
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)

test_that("footprintPerc works", {
  expect_no_error(footprintPerc(NomeData))
})

#test that footprintPerc returns the correct output with default example data
NomeData <- createExampleData(samples=c("WT_1"),group=c("WT"),nROI=10,randomMeth=FALSE)
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)
footPerc <- footprintPerc(NomeData)
footPerc$open_WT_1 <- round(footPerc$open_WT_1)

output <- tibble(ROI=c("ROI10","ROI9","ROI8","ROI7","ROI6","ROI5","ROI4","ROI3","ROI1","ROI2"),
                 ROIgroup=rep("motif1",10),
                 tf_WT_1=rep(25,10),
                 open_WT_1=rep(55,10),
                 upNuc_WT_1=rep(0,10),
                 Nuc_WT_1=rep(20,10),
                 downNuc_WT_1=rep(0,10)
)

test_that("footprintPerc returns the correct output", {
    expect_equal(footPerc,output)
   })
