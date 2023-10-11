NomeData <- createExampleData()
NomeData <- footprintCalc(NomeData)
NomeData <- footprintQuant(NomeData)
res <- diNOMeTest(NomeData)
fp <- footprintPerc(NomeData)

test_that("compareFootprints works", {
  expect_no_error(compareFootprints(fp,res,plotcols="black"))
})

