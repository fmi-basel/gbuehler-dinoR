library(tibble)

#test that it works with random example data
NomeData <- createExampleData()
test_that("metaPlots works", {
  expect_no_error(metaPlots(NomeData,span=0.08))
})

#test that the output is as expected when using default example data
NomeData <- createExampleData(samples=c("WT_1"),group=c("WT"),nROI=10,randomMeth=FALSE)
plotData <- metaPlots(NomeData,span=0.2)
plotData$protection <- round(plotData$protection)
plotData$loess <- round(plotData$loess)

output <- tibble(position = -150:150,type=rep("motif1",301),
                 protection=rep(46,301),loess=rep(46,301),
                  sample=rep("WT_1",301))

test_that("metaPlots returns the correct output", {
   expect_equal(plotData,output)
 })

