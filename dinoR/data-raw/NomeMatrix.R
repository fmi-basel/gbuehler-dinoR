## code to prepare `NomeMatrix` dataset goes here
library(tidyverse)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(nomeR)

load("/tungstenfs/scratch/gbuehler/michi/Projects/Adnp/NOMEseq/sgRNA_NOME_regions.motif.amplicons_batch1_corrected_Ctcf_enr.RData")
#combine Lys and PxVxL  and Adnp types
which$type <- ifelse(which$type=="PxVxL","Adnp",which$type)
which$type <- ifelse(which$type=="Lys","Adnp",which$type)

names(which) <- paste(which$type,seqnames(which),start(which),end(which),sep="_")

#get the amplicon sequences
seq1 <- getSeq(BSgenome.Mmusculus.UCSC.mm10,which)
#get the box subsets
seq_box1 <- subseq(seq1,261,286)
seq_box2 <- subseq(seq1,303,319)
seq_box3 <- subseq(seq1,336,361)
seq_boxes <- list(seq_box1,seq_box2,seq_box3)

#check if > 0 GCH
nmatch_per_seq_boxes <- matrix(ncol=3,nrow=length(seq1))
for (i in seq_along(seq_boxes)){
  mindex1 <- vmatchPattern("GCA", seq_boxes[[i]])
  mindex2 <- vmatchPattern("GCT", seq_boxes[[i]])
  mindex3 <- vmatchPattern("GCC", seq_boxes[[i]])
  nmatch_per_seq <- cbind(elementNROWS(mindex1),elementNROWS(mindex2),elementNROWS(mindex3))
  nmatch_per_seq_boxes[,i] <- apply(nmatch_per_seq,1,sum)
}
colnames(nmatch_per_seq_boxes) <- c("box1","box2","box3")
rownames(nmatch_per_seq_boxes) <- names(which)

#select amplicons with > 0 GCH is all 3 box GRanges files
GCHamps <- rownames(nmatch_per_seq_boxes)[apply(nmatch_per_seq_boxes,1,min)>0]

which <- which[names(which) %in% GCHamps]

#-------------------generate NomeMatrix------------------------------------------
#----------------------------------------------------------------------------
libID <- "merged"

bamdir <- paste("/tungstenfs/scratch/gbuehler/michi/Projects/Adnp/NOMEseq/bam",libID,sep=".")
bamFiles <- list.files(bamdir,pattern="_dedup.bam$",full.names = TRUE)
bamFiles <- grep("PxVxL",bamFiles,value = TRUE,invert=TRUE)
bamFiles <- grep("Ctcf",bamFiles,value = TRUE,invert=TRUE)
bamNames <- gsub(paste(bamdir,"/",sep=""),"",bamFiles)
bamNames <- gsub("_dedup.bam","",bamNames)

NomeMatrix <- get_data_matrix_from_bams(
  bamfiles= bamFiles,
  samplenames=bamNames,
  regions=which,
  genome="/tungstenfs/scratch/gbuehler/michi/Annotations/BISCUIT/mm10_withoutAltChr.fa",
  whichContext = c("GCH"),
  collapseBySample = TRUE,
  remove_nonunique = FALSE,
  clip_until_nbg = 0L,
  noclip_protect_frac_above = 0.9,
  max_bisC_meth = 0.1,
  min_bisC_size = 10,
  mapqMin = 30,
  mapqMax = 255,
  max_read_size = 1000L,
  ncores = 10L
)

NomeMatrix <- NomeMatrix[,c(1,6,10,14)]
NomeMatrix <- NomeMatrix[NomeMatrix$nFragsAnalyzed > 30 & NomeMatrix$nFragsAnalyzed < 80,]

usethis::use_data(NomeMatrix, overwrite = TRUE)
