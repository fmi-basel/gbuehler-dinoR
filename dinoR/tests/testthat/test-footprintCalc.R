NomeData <- createExampleData()

test_that("footprintCalc works", {
  expect_no_error(footprintCalc(NomeData))
})
