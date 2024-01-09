
NomeData <- createExampleData()
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)
fp <- footprintPerc(NomeData)

test_that("fpPercHeatmap works", {
  expect_no_error(fpPercHeatmap(fp))
})
