
counts <- tibble(sample = c(rep("WT",10),rep("KO",10)), ROI=rep(paste0("ROI",1:10),2),
                 TF=floor(runif(20,min=10,max=100)),open=floor(runif(20,min=10,max=100)),
                 upNuc=floor(runif(20,min=10,max=100)),Nuc=floor(runif(20,min=10,max=100)),
                 downNuc=floor(runif(20,min=10,max=100)))
test_that("footprintPerc works", {
expect_no_error(footprintPerc(counts))
})
