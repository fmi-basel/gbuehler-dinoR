#' @title diNOMeTest
#'
#' @description Tests for differential fragment counts for each NOMe
#'  footprint pattern compared to all patterns in 2 conditions.
#'
#' @details Uses edgeR's quasi-likelihood methods to conveniently test for differential proportions of each
#'  one of 5 distinct footprints between 2 control and 2 treatment samples.
#'
#' @param footprint_counts A tibble containing the sample names, ROI names,
#' and number of fragments in each NOMe footprint pattern category. For example
#' the output of the stateQuant function.
#' @param WTsamples The control sample names (two replicates) as they appear in footprint_quantifications.
#' @param KOsamples The treatment sample names (two replicates) as they appear in footprint_quantifications.
#' @param FDR The FDR cutoff for a ROI - footprint combination to be called regulated in the output.
#' @param FC The fold change cutoff for a ROI - footprint combination to be called regulated in the output.
#'
#' @return A tibble with the results of differential fragment count testing for each ROI-footprint combination.
#'
#' @examples
#' library(tibble)
#' NomeMatrix <- tibble(SampleName = c(rep("WT_1",5),
#' rep("WT_2",5),rep("KO_1",5),rep("KO_2",5)),
#' names=rep(paste0("ROI",1:5),4),nFragsAnalyzed=rep(20,20),
#' GCH_DataMatrix=rep(list(matrix(sample(c(0,1),size=150*20,
#' replace=TRUE),ncol=150,nrow=20)),20))
#' patternQuant <- footprintQuant(NomeMatrix)
#' diNOMeTest(patternQuant,WTsamples = c("WT_1","WT_2"),
#' KOsamples = c("KO_1","KO_2"))
#'
#' @importFrom stats model.matrix
#' @importFrom stats complete.cases
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom edgeR DGEList
#' @importFrom edgeR calcNormFactors
#' @importFrom edgeR estimateDisp
#' @importFrom edgeR glmQLFit
#' @importFrom edgeR glmQLFTest
#' @importFrom edgeR topTags
#'
#'
#' @export
diNOMeTest <- function(footprint_counts,WTsamples = c("WT_1","WT_2"),
                       KOsamples = c("KO_1","KO_2"),
                       FDR=0.05,FC=2){
  #model
  annots <- data.frame(cond= c(rep("WT",12),rep("KO",12)),
                     rep= c(c(rep("r1",6),rep("r2",6)),rep("r3",6),rep("r4",6)),
                     type= c(rep(c("all","TF","open","Nuc","upNuc","downNuc"),4)))

  dm <- model.matrix(~ rep,data=annots)
  dm <- cbind(dm,openWT = as.numeric(annots$type =="open" & annots$cond =="WT"))
  dm <- cbind(dm,openKO = as.numeric(annots$type =="open" & annots$cond =="KO"))
  dm <- cbind(dm,TFWT = as.numeric(annots$type =="TF" & annots$cond =="WT"))
  dm <- cbind(dm,TFKO = as.numeric(annots$type =="TF" & annots$cond =="KO"))
  dm <- cbind(dm,nucWT = as.numeric(annots$type =="Nuc" & annots$cond =="WT"))
  dm <- cbind(dm,nucKO = as.numeric(annots$type =="Nuc" & annots$cond =="KO"))
  dm <- cbind(dm,dnucWT = as.numeric(annots$type =="downNuc" & annots$cond =="WT"))
  dm <- cbind(dm,dnucKO = as.numeric(annots$type =="downNuc" & annots$cond =="KO"))
  dm <- cbind(dm,unucWT = as.numeric(annots$type =="upNuc" & annots$cond =="WT"))
  dm <- cbind(dm,unucKO = as.numeric(annots$type =="upNuc" & annots$cond =="KO"))

  #define contrasts
  contrast_open_vs_all <- c(0,0,0,0,-1,1,0,0,0,0,0,0,0,0)
  contrast_TF_vs_all <- c(0,0,0,0,0,0,-1,1,0,0,0,0,0,0)
  contrast_nuc_vs_all <- c(0,0,0,0,0,0,0,0,-1,1,0,0,0,0)
  contrast_dnuc_vs_all <- c(0,0,0,0,0,0,0,0,0,0,-1,1,0,0)
  contrast_unuc_vs_all <- c(0,0,0,0,0,0,0,0,0,0,0,0,-1,1)

  contrasts <- list(contrast_open_vs_all,contrast_TF_vs_all,contrast_nuc_vs_all,contrast_dnuc_vs_all,contrast_unuc_vs_all)
  names(contrasts) <- c("open_vs_all","TF_vs_all","Nuc_vs_all","downNuc_vs_all","upNuc_vs_all")

  #model using EdgeR and design matrix from above
  #re arrange pattern quantification matrix
  footprint_counts$all <- apply(footprint_counts[,3:7],1,sum)
  footprint_counts <- pivot_wider(footprint_counts,id_cols="ROI",
                              names_from = .data$sample,values_from=c(.data$all,.data$TF,.data$open,
                                                                      .data$Nuc,.data$upNuc,.data$downNuc))

  #only keep ROIs where there were data in all samples
  footprint_counts <- footprint_counts[complete.cases(footprint_counts)==TRUE,]
  footprint_counts2 <- as.matrix(footprint_counts[,-1])
  row.names(footprint_counts2) <- footprint_counts$ROI

  #add sample names to annots using the specifications which samples are WT and which samples are KO
  annots$sample <- c( paste(c("all","TF","open","Nuc","upNuc","downNuc"),WTsamples[1],sep="_"),
  paste(c("all","TF","open","Nuc","upNuc","downNuc"),WTsamples[2],sep="_"),
  paste(c("all","TF","open","Nuc","upNuc","downNuc"),KOsamples[1],sep="_"),
  paste(c("all","TF","open","Nuc","upNuc","downNuc"),KOsamples[2],sep="_"))

  #order pattern quantifiction matrix by annots sample names
  footprint_counts2 <- footprint_counts2[,annots$sample]

  #calculate combined library sizes based on total patterned reads in sample
  libSizes1 <- colSums(footprint_counts2)
  libSizes2 <- c(rep(sum(libSizes1[grep(WTsamples[1],names(libSizes1))]),6),
                 rep(sum(libSizes1[grep(WTsamples[2],names(libSizes1))]),6),
                 rep(sum(libSizes1[grep(KOsamples[1],names(libSizes1))]),6),
                 rep(sum(libSizes1[grep(KOsamples[2],names(libSizes1))]),6))
  names(libSizes2) <- names(libSizes1)

  #generate PatternQuant matrix with read counts per sample
  footprint_counts_sums <- as.matrix(data.frame(WT_1 = apply(footprint_counts2[,1:6],1,sum),
                                           WT_2 = apply(footprint_counts2[,7:12],1,sum),
                                           KO_1 = apply(footprint_counts2[,13:18],1,sum),
                                           KO_2 = apply(footprint_counts2[,19:24],1,sum)))

  #calculate lib size and norm factors on sample wise matrix
  ySums <- DGEList(counts=footprint_counts_sums)
  ySums <- calcNormFactors(ySums,method="TMM")
  SummedNormFactors <- c(rep(ySums$samples$norm.factors[1],6),rep(ySums$samples$norm.factors[2],6),
                         rep(ySums$samples$norm.factors[3],6),rep(ySums$samples$norm.factors[4],6))

  #generate DGE list object for each sample/type combination
  y <- DGEList(counts=footprint_counts2,lib.size=libSizes2) #add lib sizes based on complete read number per sample

  #add norm factors calculated on Summed values
  y$samples$norm.factors <- SummedNormFactors
  #dispersion
  y <- estimateDisp(y,dm)
  #model
  fit <- glmQLFit(y, dm)

  # use the contrasts defined above to get p-values and fold-changes
  res <- list()
  for (i in seq_along(contrasts)){
    qlf <- glmQLFTest(fit, contrast=contrasts[[i]])
    res[[i]] <- data.frame(topTags(qlf,n=nrow(footprint_counts2),adjust.method = "BH"))
    res[[i]]$contrasts <- names(contrasts)[i]
    res[[i]]$ROI <- row.names(res[[i]])
  }
  res2 <- do.call("rbind",res)
  res2$logadjPval <- -log10(res2$FDR)
  res2$regulated <- ifelse(res2$FDR < FDR & res2$logFC > log2(FC),"up",
                           ifelse(res2$FDR < FDR & res2$logFC < -log2(FC),"down","no"))
  return(tibble(res2))


}
