#' @title diNOMeTest
#'
#' @description Tests for differences in fragment counts for each NOMe
#'  footprint pattern compared to total fragment counts in two conditions.
#'
#' @details Uses edgeR's quasi-likelihood methods to conveniently test for differential proportions of each
#'  one of 5 (or 3, if nucleosome footprints are combined) distinct footprints between at least two control and at least two treatment samples.
#'
#' @param footprint_counts A Summarized Experiment containing the sample names (colData), ROI names (rowData),
#' and number of fragments in each NOMe footprint pattern category (assays). For example
#' the output of the (\code{footprintQuant}) function.
#' @param WTsamples The control sample names as they appear in (\code{colData(footprint_counts)$samples}).
#' @param KOsamples The treatment sample names as they appear in (\code{colData(footprint_counts)$samples}).
#' @param minreads The minimum number of fragments to which a footprint could be assigned
#'  a ROI must have in all samples. All other ROIs are filtered out before the differential NOMe analysis.
#' @param meanreads The minimum number of fragments to which a footprint could be assigned
#'  a ROI must on average across all samples. All other ROIs are filtered out before the differential NOMe analysis.
#' @param prior.count The pseudocount used for (\code{edgeR::glmQLFit}).
#' @param FDR The FDR cutoff for a ROI - footprint combination to be called regulated in the output.
#' @param FC The fold change cutoff for a ROI - footprint combination to be called regulated in the output.
#' @param combineNucCounts If TRUE, the upNuc, downNuc, and Nuc fragment counts will be combined into the Nuc category.
#'
#' @return A tibble with the results of differential fragment count testing for each ROI-footprint combination.
#'
#' @examples
#' NomeData <- createExampleData()
#' NomeData <- footprintCalc(NomeData)
#' footprint_counts <- footprintQuant(NomeData)
#' diNOMeTest(footprint_counts,WTsamples = c("WT_1","WT_2"),
#' KOsamples = c("KO_1","KO_2"))
#'
#' @importFrom stats model.matrix complete.cases
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom SummarizedExperiment assays rowData
#'
#' @export
diNOMeTest <- function(footprint_counts,WTsamples=c("WT_1","WT_2"),
                       KOsamples=c("KO_1","KO_2"), minreads=1,meanreads=1, prior.count=3,
                       FDR=0.05,FC=2,combineNucCounts=FALSE){

  #keep only ROIs where all samples have more than 0 total counts
  footprint_counts <- footprint_counts[which(apply(assays(footprint_counts)[["all"]],
                                                   1,min,na.rm=TRUE) >= minreads),]
  footprint_counts <- footprint_counts[which(apply(assays(footprint_counts)[["all"]],
                                                   1,mean,na.rm=TRUE) >= meanreads),]

  #determine number of replicates
  nrepWT <- length(WTsamples)
  nrepKO <- length(KOsamples)
  nreps <- nrepWT + nrepKO
  allSamples <- c(WTsamples,KOsamples)

  #check that tehre are at least 2 replicates for each condition
  if (nrepWT < 2 | nrepKO < 2){
    stop("Please provide at least two replicate samples for each condition!")
  }

  #model
  if (combineNucCounts == TRUE){
    annots <- data.frame(cond= c(rep("WT",4*nrepWT),rep("KO",4*nrepKO)),
                         rep= rep(paste0("r",seq(1,nreps)),each=4),
                         type= c(rep(c("all","tf","open","Nuc"),nreps)))
  } else {
  annots <- data.frame(cond= c(rep("WT",6*nrepWT),rep("KO",6*nrepKO)),
                     rep= rep(paste0("r",seq(1,nreps)),each=6),
                     type= c(rep(c("all","tf","open","Nuc","upNuc","downNuc"),nreps)))
  }

  dm <- model.matrix(~ rep,data=annots)
  dm <- cbind(dm,openWT = as.numeric(annots$type =="open" & annots$cond =="WT"))
  dm <- cbind(dm,openKO = as.numeric(annots$type =="open" & annots$cond =="KO"))
  dm <- cbind(dm,tfWT = as.numeric(annots$type =="tf" & annots$cond =="WT"))
  dm <- cbind(dm,tfKO = as.numeric(annots$type =="tf" & annots$cond =="KO"))
  dm <- cbind(dm,nucWT = as.numeric(annots$type =="Nuc" & annots$cond =="WT"))
  dm <- cbind(dm,nucKO = as.numeric(annots$type =="Nuc" & annots$cond =="KO"))
  if (combineNucCounts == FALSE){
  dm <- cbind(dm,dnucWT = as.numeric(annots$type =="downNuc" & annots$cond =="WT"))
  dm <- cbind(dm,dnucKO = as.numeric(annots$type =="downNuc" & annots$cond =="KO"))
  dm <- cbind(dm,unucWT = as.numeric(annots$type =="upNuc" & annots$cond =="WT"))
  dm <- cbind(dm,unucKO = as.numeric(annots$type =="upNuc" & annots$cond =="KO"))
  }

  #define contrasts
  if (combineNucCounts == FALSE){
    contrast_open_vs_all <- c(rep(0,nreps),-1,1,0,0,0,0,0,0,0,0)
    contrast_tf_vs_all <- c(rep(0,nreps),0,0,-1,1,0,0,0,0,0,0)
    contrast_nuc_vs_all <- c(rep(0,nreps),0,0,0,0,-1,1,0,0,0,0)
    contrast_dnuc_vs_all <- c(rep(0,nreps),0,0,0,0,0,0,-1,1,0,0)
    contrast_unuc_vs_all <- c(rep(0,nreps),0,0,0,0,0,0,0,0,-1,1)

    contrasts <- list(contrast_open_vs_all,contrast_tf_vs_all,contrast_nuc_vs_all,
                    contrast_dnuc_vs_all,contrast_unuc_vs_all)
    names(contrasts) <- c("open_vs_all","tf_vs_all","Nuc_vs_all","downNuc_vs_all","upNuc_vs_all")
  } else {
    contrast_open_vs_all <- c(rep(0,nreps),-1,1,0,0,0,0)
    contrast_tf_vs_all <- c(rep(0,nreps),0,0,-1,1,0,0)
    contrast_nuc_vs_all <- c(rep(0,nreps),0,0,0,0,-1,1)

    contrasts <- list(contrast_open_vs_all,contrast_tf_vs_all,contrast_nuc_vs_all)
    names(contrasts) <- c("open_vs_all","tf_vs_all","Nuc_vs_all")
  }
  #model using EdgeR and design matrix from above
  #re arrange pattern quantification matrix from summarized experiment
  patterns <- unique(annots$type)
  assay_list <- list()
  for (p in seq_along(patterns)){
    assay_list[[p]] <- assays(footprint_counts)[[patterns[p]]]
    colnames(assay_list[[p]]) <- paste(patterns[p], colnames(assay_list[[p]]),sep="_")
  }
  footprint_counts1 <- data.frame(ROI=rownames(rowData(footprint_counts)),
                                  do.call("cbind",assay_list))

  #only keep ROIs where there were data in all samples
  footprint_counts1 <- footprint_counts1[complete.cases(footprint_counts1)==TRUE,]
  footprint_counts2 <- as.matrix(footprint_counts1[,-1])
  row.names(footprint_counts2) <- footprint_counts1$ROI

  #add sample names to annots using the specifications which samples are WT and which samples are KO
  sample_names <- character(0)
   for (i in seq_along(allSamples)){
    sample_names <- c(sample_names,paste(patterns,allSamples[i],sep="_"))
  }
  annots$sample <- sample_names

  #order pattern quantifiction matrix by annots sample names
  footprint_counts2 <- footprint_counts2[,annots$sample]

  #generate PatternQuant matrix with read counts per sample
  footprint_counts_sums <- matrix(ncol=length(allSamples),nrow=nrow(footprint_counts2))
  for (i in seq_along(allSamples)){
    footprint_counts_sums[,i] <-  apply(footprint_counts2[,grep(allSamples[i],
                                                                colnames(footprint_counts2))],1,sum)
  }
  colnames(footprint_counts_sums) <- allSamples
  rownames(footprint_counts_sums) <- rownames(footprint_counts2)

  #calculate lib size and norm factors on sample wise matrix
  ySums <- DGEList(counts=footprint_counts_sums)
  ySums <- calcNormFactors(ySums,method="TMM")
  libSizes <- rep(ySums$samples$lib.size,each=length(patterns))
  SummedNormFactors <- rep(ySums$samples$norm.factors,each=length(patterns))

  #generate DGE list object for each sample/type combination...
  #...adding lib sizes based on complete read number per sample
  y <- DGEList(counts=footprint_counts2,lib.size=libSizes)

  #...adding norm factors calculated on Summed values
  y$samples$norm.factors <- SummedNormFactors

  #dispersion
  y <- estimateDisp(y,dm)
  #model
  fit <- glmQLFit(y, dm,prior.count=prior.count)

  # use the contrasts defined above to get p-values and fold-changes
  res <- list()
  for (i in seq_along(contrasts)){
    qlf <- glmQLFTest(fit, contrast=contrasts[[i]])
    res[[i]] <- data.frame(topTags(qlf,n=nrow(footprint_counts2),adjust.method = "BH"))
    res[[i]]$contrasts <- names(contrasts)[i]
    res[[i]]$ROI <- row.names(res[[i]])
    res[[i]] <- cbind(res[[i]],rowData(footprint_counts))
  }
  res2 <- do.call("rbind",res)
  res2$logadjPval <- -log10(res2$FDR)
  res2$regulated <- ifelse(res2$FDR < FDR & res2$logFC > log2(FC),"up",
                           ifelse(res2$FDR < FDR & res2$logFC < -log2(FC),"down","no"))
  return(tibble(res2))
}
