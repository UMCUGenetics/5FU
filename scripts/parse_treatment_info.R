library("readxl")
library(stringr)
library(tidyr)
library(plyr)
library(dplyr)
library(devtools)
#install_github("UMCUGenetics/MutationalPatterns", ref = "develop")
library("MutationalPatterns")
#biocLite("dendextend")
library("GenomeInfoDb")
#biocLite("xlsx")
library(BSgenome)
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(reshape2)
library(ggplot2) 
library(MutationalPatterns)
library(GenomicRanges)
library(VariantAnnotation)
library(ggdendro)
library(cowplot)
library(dendextend)
library(stringr)
library(xlsx)
library("readxl")
library(dplyr)
library(foreach)
library(caret)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(dunn.test)
library(FSA)

#bugfix
#remove DelayedArray()
#install_github("Bioconductor-mirror/DelayedArray", ref="22cf715874e9c9378b9ff3691d7cb72fc82bcbf6")
library(DelayedArray)

# sessionInfo()
# packageVersion("NMF")

#check R version
# version
# packageStatus()
# packageVersion("DelayedArray")

options(stringsAsFactors = F)

# External functions

source('~/Documents/HMF_data/scripts/R_functions.R')
source('~/Documents/HMF_data/scripts/plotOverviewAllSignatures.R')
source('~/Documents/HMF_data/scripts/plotContributionPerSignature.R')
source('~/Documents/HMF_data/scripts/plotBoxPlotTreatments.R')
source('~/Documents/HMF_data/scripts/loadMutationMatrix.R')
source('~/Documents/HMF_data/scripts/loadTGCounts.R')

###########################################################################################################################
########################## DATA LOADING ###################################################################################

# Default colours for mutation spectrum plotting
COLORS6 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#D4D2D2", "#ADCC54", "#F0D0CE")

COLORS7 = c(
  "#2EBAED", "#000000", "#DE1C14",
  "#E98C7B", "#D4D2D2", "#ADCC54",
  "#F0D0CE")

SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
SUBSTITUTIONS_96 = rep(SUBSTITUTIONS, each=16)
SUBSTITUTIONS_192 = rep(SUBSTITUTIONS, each=32)

C_TRIPLETS = c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS = c(
  "ATA", "ATC", "ATG", "ATT",
  "CTA", "CTC", "CTG", "CTT",
  "GTA", "GTC", "GTG", "GTT",
  "TTA", "TTC", "TTG", "TTT")

TRIPLETS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

STRAND = rep(c("U","T"), 96)
DNA_BASES = c("A", "C", "G", "T")

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# cancer_col=c("Lung"="#7FC97F",
#              "Esophagus"="#BEAED4",
#              "Breast"="#FDC086",
#              "Skin" ="#FFFF99",
#              "Ovary" ="#386CB0",
#              "Prostate" ="#F0027F",
#              "Biliary" ="#BF5B17",
#              "Testis" = "#666666",
#              "Colon/Rectum" ="#1B9E77",
#              "Bone/Soft tissue"= "#D95F02",
#              "Liver" ="#7570B3",
#              "NET" ="#E7298A",
#              "Kidney" ="#66A61E",
#              "Pancreas" ="#E6AB02",
#              "Head and neck"= "#A6761D",
#              "Unknown" ="#A6CEE3",
#              "Thyroid" ="#1F78B4",
#              "Urinary tract" ="#B2DF8A",
#              "Penile"= "#FB9A99",
#              "Uterus" ="#E31A1C",
#              "Vulva" ="#FDBF6F",
#              "Mesothelioma"= "#FF7F00",
#              " "= "#6A3D9A",
#              "Stomach" ="#B15928",
#              "Lymphoid"= "#FBB4AE",
#              "Nervous system"= "#B3CDE3",
#              "CUP"= "#CCEBC5",
#              "Small intestine"= "#DECBE4",
#              "Adrenal"= "#FED9A6",
#              "Double primary"="#FFFFCC",
#              "Thymus"= "#E5D8BD",
#              "Myeloid"= "#FDDAEC",
#              "Eye" ="#F2F2F2")

colorpalette2=c("SBS1"="#7FC97F",
                  "SBS2"="#BEAED4",
                  "SBS3"="#FDC086",
                  "SBS4"="#FFFF99",
                  "SBS5"="#386CB0",
                  "SBS6"="#F0027F",
                  "SBS7a"="#BF5B17",
                  "SBS7b"="#666666",
                  "SBS7c"="#1B9E77",
                  "SBS7d"="#D95F02",
                  "SBS8"="#7570B3",
                  "SBS9"="#E7298A",
                  "SBS10a"="#66A61E",
                  "SBS10b"="#E6AB02",
                  "SBS11"="#A6761D",
                  "SBS12"="#A6CEE3",
                  "SBS13"="#1F78B4",
                  "SBS14"="#B2DF8A",
                  "SBS15"="#33A02C",
                  "SBS16"="#FB9A99",
                  "SBS17a"="#E31A1C",
                  "SBS17b"="#FDBF6F",
                  "SBS18"="#FF7F00",
                  "SBS19"="#CAB2D6",
                  "SBS20"="#6A3D9A",
                  "SBS21"="#B15928",
                  "SBS22"="#FBB4AE",
                  "SBS23"="#B3CDE3",
                  "SBS24"="#CCEBC5",
                  "SBS25"="#DECBE4",
                  "SBS26"="#FED9A6",
                  "SBS27"="#FFFFCC",
                  "SBS28"="#E5D8BD",
                  "SBS29"="#FDDAEC",
                  "SBS30"="#F2F2F2",
                  "SBS31"="#B3E2CD",
                  "SBS32"="#FDCDAC",
                  "SBS33"="#CBD5E8",
                  "SBS34"="#F4CAE4",
                  "SBS35"="#E6F5C9",
                  "SBS36"="#FFF2AE",
                  "SBS37"="#F1E2CC",
                  "SBS38"="#CCCCCC",
                  "SBS39"="#E41A1C",
                  "SBS40"="#377EB8",
                  "SBS41"="#4DAF4A",
                  "SBS42"="#984EA3",
                  "SBS43"="#FFFF33",
                  "SBS44"="#A65628",
                  "SBS45"="#F781BF",
                  "SBS46"="#999999",
                  "SBS47"="#66C2A5",
                  "SBS48"="#FC8D62",
                  "SBS49"="#8DA0CB",
                  "SBS50"="#E78AC3",
                  "SBS51"="#A6D854",
                  "SBS52"="#FFD92F",
                  "SBS53"="#E5C494",
                  "SBS54"="#B3B3B3",
                  "SBS55"="#8DD3C7",
                  "SBS56"="#FFFFB3",
                  "SBS57"="#BEBADA",
                  "SBS58"="#FB8072",
                  "SBS59"="#80B1D3",
                  "SBS60"="#FDB462")

colors_treatment=c("<NA>"="#FDB462",
                   "1"="#7FC97F",
                   "2"="#BEAED4",
                   "3"="#FDC086",
                   "4"="#FFFF99",
                   "5"="#386CB0",
                   "6"="#F0027F",
                   "7"="#D95F02",
                   "8"="#7570B3",
                   "9"="#E7298A",
                   "10"="#66A61E",
                   "10"="#E6AB02",
                   "11"="#A6761D",
                   "12"="#A6CEE3",
                   "13"="#1F78B4",
                   "14"="#B2DF8A",
                   "15"="#33A02C")

clinical <- data.frame(read_excel("~/Documents/HMF_data/treatment_info/clinical.xlsx", sheet='clinical',na="NULL"))
clinical_DR10 <- subset(clinical, DR.10.update!='NA')
names(clinical_DR10)
clinical_DR10$preTreatments_dub=clinical_DR10$preTreatments
clinical_DR10$preTreatmentsType_dub=clinical_DR10$preTreatmentsType
max_fields=max(sapply(strsplit(as.character(clinical_DR10$preTreatments_dub),'/'),length))
clinical_DR10=separate(data=clinical_DR10, col=preTreatments_dub, into=paste0("preTreatment_",1:as.numeric(max_fields)), sep="/")
clinical_DR10=separate(data=clinical_DR10, col=preTreatmentsType_dub, into=paste0("preTreatmentsType_",1:as.numeric(max_fields)), sep="/")
clinical_DR10_pretreats=clinical_DR10[paste0("preTreatment_",1:as.numeric(max_fields))]
clinical_DR10_pretreats_named=cbind(clinical_DR10$sampleId,clinical_DR10_pretreats)
clinical_DR10_preTreatmentsType=clinical_DR10[paste0("preTreatmentsType_",1:as.numeric(max_fields))]

uq_elem=c()
for(i in 1:ncol(clinical_DR10_pretreats))
{
  uq_elem=c(unique(clinical_DR10_pretreats[,i]), uq_elem)
  uq_elem=unique(uq_elem)
  uq_elem=uq_elem[!is.na(uq_elem)]
}

treatment_type=c()
for(i in 1:ncol(clinical_DR10_pretreats))
{
    treatment_type_sub=paste(clinical_DR10_pretreats[,i],clinical_DR10_preTreatmentsType[,i],sep="_")
    treatment_type=c(unique(treatment_type_sub), treatment_type)
    treatment_type=unique(treatment_type)
}

df_final = setNames(data.frame(matrix(ncol = length(uq_elem)+1, nrow = 0)),c("sampleID",uq_elem))
for(i in 1:nrow(clinical_DR10_pretreats)) {
  row=clinical_DR10_pretreats_named[i,]
  row=as.character(unlist(as.vector(lapply(row, as.character))))
  df1 = setNames(data.frame(matrix(ncol = length(uq_elem)+1, nrow = 0)),c("sampleID",uq_elem))
  names(df1)
  
  name=row[1]
  if(sum(nchar(row[!is.na(row)][-1]))>0){
    df2=as.data.frame(t(unstack(as.data.frame(table(row[row %in% uq_elem])), form=Freq~Var1)))
    df3=dplyr::full_join(df1, df2)
    df3$sampleID=name
  }else{
    df3=setNames(data.frame(matrix(ncol = length(uq_elem)+1, nrow = 1)),c("sampleID",uq_elem))
    df3$sampleID=name
  }
  df_final=rbind(df_final,df3)
  rm(df1,df2,df3)
}

non_treated <- vector()
for (i in 1:nrow(df_final)){
  if (all(is.na(df_final[i,2:191]))){
    non_treated <- c(non_treated, i)
  }
}
df_final <- df_final[-non_treated,]

# View(head(df_final))
metadata <- data.frame(read.csv("~/Documents/HMF_data/Metadata_DR10-update/DR-010-update_metadata_180815.tsv",header = T,sep = "\t"))
df_final$CancerType <- unlist(lapply(df_final$sampleID, function(x){
  x <- factor(x, levels=unique(metadata$sampleId))
  if (is.na(x)){
    return(NA)
  } else {
    rowNum <- which(metadata$sampleId==x)
    if (length(rowNum)==0){
      return(NA)
    } else {
      type <- as.character(metadata$primaryTumorLocation[rowNum])
      if (type == "Bone/soft tissue"){
        type <- "Bone/Soft tissue"
      } else if (type == "Head and Neck"){
        type <- "Head and neck"
      }
      return(type)
    }
  }
}))

### analyse all data
##HMF_data
clonal <- "SUBCLONAL"
mut_matrix_all <- loadMutationMatrix(path="~/Documents/HMF_data/Metadata_DR10-update/",
                                     pattern="96.txt",
                                     clonality = clonal)

# mut_mat = mut_matrix_all + 0.0001

TG_counts <- loadTGCounts(path="~/Documents/HMF_data/Metadata_DR10-update/",
                          pattern="96.txt",
                          clonality = clonal)

cancer_signatures_new = read.csv("~/Documents/HMF_data/matrices/sigProfiler_SBS_signatures.csv", sep = ";", header = TRUE)
cancer_signatures_new = cancer_signatures_new[order(cancer_signatures_new[,1]),]
cancer_signatures_new = as.matrix(cancer_signatures_new[,4:68])

fit_res <- fit_to_signatures(mut_mat, cancer_signatures_new)
contribution_cosm = fit_res$contribution

cancer_col <- unique(col_vector)[1:length(unique(metadata$primaryTumorLocation))]
names(cancer_col) <- unique(metadata$primaryTumorLocation)

###########################################################################################################################
####################################### MAKE PLOTS #######################################################################

system.time({
plotOverviewAllSignatures(fc = contribution_cosm, df_treatment = df_final, 
                          color_sign = colorpalette2, color_treat = colors_treatment, color_cancer = cancer_col, 
                          signature = "SBS17b", treatment = "Fluorouracil", 
                          clonal="CLONAL", TG = TG_counts)
  
plotOverviewAllSignatures(fc = contribution_cosm, df_treatment = df_final, 
                          color_sign = colorpalette2, color_treat = colors_treatment, color_cancer = cancer_col, 
                          signature = "SBS17b", treatment = "Fluorouracil", abs=F, clonal="CLONAL")

plotContributionPerSignature(fc = contribution_cosm, df_treatment = df_final, 
                             color_sign = colorpalette2, color_treat = colors_treatment, 
                             signature = "SBS17b", treatment = "Fluorouracil", clonal="CLONAL")

plotBoxPlotTreatments(fc = contribution_cosm, df_treatment = df_final, 
                      color_sign = colorpallete2, color_treat=colors_treatment, 
                      signature = "SBS17b", treatment="Fluorouracil", 
                      showCancerType = T, showTreatment = T, showClonality =T,TG=TG_counts, clonal = clonal)

plotBoxPlotTreatments(fc = contribution_cosm, df_treatment = df_final,
                      color_sign = colorpallete2, color_treat=colors_treatment, 
                      signature = "SBS17b", treatment="Fluorouracil", 
                      showCancerType = T, showTreatment = F)

plotBoxPlotTreatments(fc = contribution_cosm, df_treatment = df_final, 
                      color_sign = colorpallete2, color_treat=colors_treatment, 
                      signature = "SBS17b", treatment="Fluorouracil", 
                      showCancerType = F, showTreatment = T)

plotBoxPlotTreatments(fc = contribution_cosm, df_treatment = df_final, 
                      color_sign = colorpallete2, color_treat=colors_treatment, 
                      signature = "SBS17b", treatment="Fluorouracil", 
                      showCancerType = F, showTreatment = F, TG = TG_counts)
})
