
library(tibble)
NomeMatrix <- tibble(SampleName = c(rep("WT_1",5),
rep("WT_2",5),rep("KO_1",5),rep("KO_2",5)),
names=rep(paste0("ROI",1:5),4),nFragsAnalyzed=rep(20,20),
GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
replace=TRUE),ncol=150,nrow=20)),20))
footprint_counts <- footprintQuant(NomeMatrix)

test_that("diNOMeTest works", {
expect_no_error(diNOMeTest(footprint_counts,WTsamples = c("WT_1","WT_2"),
KOsamples = c("KO_1","KO_2")))
})
