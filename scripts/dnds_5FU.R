library(foreach)
library(dndscv)


breast_5FU <- list.files("~/surfdrive/Shared/Sig17/HMF_data/dnds/5FU/breast_5FU/", recursive = T, full.names = T)
colon_5FU <- list.files("~/surfdrive/Shared/Sig17/HMF_data/dnds/5FU/colon_5FU/", recursive = T, full.names = T)

mutations_5FU <- foreach(f = c(breast_5FU,colon_5FU), .combine = 'rbind') %do% {
  df <- read.csv(f, sep= " ")
  return(df)
}

mutations_5FU <- mutations_5FU[,c("sampleid","chromosome","position", "ref","alt")]
colnames(mutations_5FU) <- c("sampleID", "chr", "pos", "ref", "mut")
mutations_5FU$chr <- substr(mutations_5FU$chr, 4, nchar(as.vector(mutations_5FU$chr)))
dndsout_5FU = dndscv(mutations_5FU)
dndsout_5FU$annotmuts
sel_cv_5FU = dndsout_5FU$sel_cv
signif_genes_5FU = sel_cv_5FU[sel_cv_5FU$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

## breast_not5FU
breast_5FU <- list.files("~/surfdrive/Shared/Sig17/HMF_data/dnds/5FU/breast_not5FU/", recursive = T, full.names = T)
colon_5FU <- list.files("~/surfdrive/Shared/Sig17/HMF_data/dnds/5FU/colon_not5FU/", recursive = T, full.names = T)
mutations_not5FU <- foreach(f = c(breast_not5FU,colon_not5FU), .combine = 'rbind') %do% {  
  df <- read.csv(f, sep= " ")
  return(df)
}

mutations_not5FU <- mutations_not5FU[,c("sampleid","chromosome","position", "ref","alt")]
colnames(mutations_not5FU) <- c("sampleID", "chr", "pos", "ref", "mut")

mutations_not5FU$chr <- substr(mutations_not5FU$chr, 4, nchar(as.vector(mutations_not5FU$chr)))

dndsout_not5FU = dndscv(mutations_not5FU)
sel_cv_not5FU = dndsout_not5FU$sel_cv
signif_genes_not5FU = sel_cv_not5FU[sel_cv_not5FU$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

## Significant genes only positive selected after 5FU treatment
genes_5FU_specific <- signif_genes_5FU$gene_name[which(is.na(match(signif_genes_5FU$gene_name, signif_genes_not5FU$gene_name)))]
length(unique(dndsout_not5FU$annotmuts[which(dndsout_not5FU$annotmut$gene=="BCL9L"),]$sampleID))
length(unique(dndsout_not5FU$annotmuts$sampleID))
length(unique(dndsout_5FU$annotmuts[which(dndsout_5FU$annotmut$gene=="BCL9L"),]$sampleID))
length(unique(dndsout_5FU$annotmuts$sampleID))


