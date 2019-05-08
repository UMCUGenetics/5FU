library(dplyr)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)


########################
########################
####COLON + BREAST + invitro
########################
########################
load("~/surfdrive/Shared/Sig17/HMF_data/cohort_analyse/somatic_clinical.RData")
dirpath="~/surfdrive/Shared/Sig17/HMF_data/plots/DR47/COLON_BREAST_INVITRO/COLON/"
dirpath="~/surfdrive/Shared/Sig17/HMF_data/plots/DR47/COLON_BREAST_INVITRO/BREAST/"
plot_name="denovo_NMF_16"


#genome_length <- sum(as.numeric(as.vector(as.data.frame(seqlengths(Hsapiens))[1:24,])))/1000000
genome_length <- 2858674662/1000000
somatic_clinical$TMB = somatic_clinical$mut_load/genome_length
somatic_clinical_subset = somatic_clinical %>% 
  filter(primaryTumorLocation == "Breast",Fluorouracil != "-2",TMB <=10) #Colon/Rectum #Breast

somatic_clinical_subset %>% group_by(Fluorouracil) %>%
  dplyr::summarize(count = n(), 
                   mean_sign = mean(NMF_H),median_sign = median(NMF_H), SD_sign = sd(NMF_H),
                   sum_sign=sum(NMF_H),max_sign=max(NMF_H),min_sign=min(NMF_H),
                   mean_mut_load = mean(mut_load), SD_mut_load = sd(mut_load),
                   median_mut_load = median(mut_load),
                   sum_mut_load=sum(mut_load))

names(somatic_clinical)


FU <- somatic_clinical_subset  %>% filter(Fluorouracil == "1") %>% pull(sampleId)
nonFU <- somatic_clinical_subset  %>% filter(Fluorouracil == "0") %>% pull(sampleId)

DATA_nonsynonymous <- read.table("~/surfdrive/Shared/Sig17/HMF_data/colon_nonsynonymous_MoA.txt", skip = 1, stringsAsFactors = FALSE, sep = " ")
DATA_nonsynonymous <- DATA_nonsynonymous[, 1:26]
head(DATA_nonsynonymous)
DATA_synonymous <- read.table("~/surfdrive/Shared//Sig17/HMF_data/colon_synonymous_MoA.txt", skip = 1, stringsAsFactors = FALSE, sep = " ")
DATA_synonymous <- DATA_synonymous[, 1:26]
nrow(DATA_nonsynonymous)

DATA_nonsynonymous <- read.table("~/surfdrive/Shared//Sig17/HMF_data/breast_nonsynonymous_MoA.txt", skip = 1, stringsAsFactors = FALSE, sep = " ")
DATA_nonsynonymous <- DATA_nonsynonymous[, 1:27]
head(DATA_nonsynonymous)
DATA_synonymous <- read.table("~/surfdrive/Shared//Sig17/HMF_data/breast_synonymous_MoA.txt", skip = 1, stringsAsFactors = FALSE, sep = " ")
DATA_synonymous <- DATA_synonymous[, 1:27]

table(DATA_nonsynonymous[which(DATA_nonsynonymous$clonality=="SUBCLONAL"),]$driver)

#fisher on drivers
#treated or not
pathos_nonFU <- nrow(DATA_nonsynonymous %>% filter(sampleid%in% nonFU, driver=="Yes",clonality=="SUBCLONAL") %>% unique())
samples_nonFU <- length(DATA_nonsynonymous %>% filter(sampleid%in% nonFU) %>% pull(sampleid) %>% unique())
pathos_FU <- nrow(DATA_nonsynonymous %>% filter(sampleid%in% FU, driver=="Yes",clonality=="SUBCLONAL") %>% unique())
samples_FU <-  length(DATA_nonsynonymous %>% filter(sampleid%in% FU) %>% pull(sampleid) %>% unique())

challenge.df = matrix(c(pathos_FU,pathos_nonFU,samples_FU,samples_nonFU), nrow = 2)
fisher.test(challenge.df, alternative ="greater")

pathos_FU_yes <- nrow(DATA_nonsynonymous %>% filter(sampleid%in% FU, driver=="Yes",clonality=="CLONAL"))
pathos_FU_no <- nrow(DATA_nonsynonymous %>% filter(sampleid%in% FU, driver!="Yes",clonality=="CLONAL"))
print(pathos_FU_yes/pathos_FU_no)

syn_5FU <- nrow(DATA_synonymous %>% filter(sampleid%in% FU, clonality=="SUBCLONAL",NMF_H>0.5))
nonsyn_5FU <- nrow(DATA_nonsynonymous %>% filter(sampleid%in% FU, clonality=="SUBCLONAL",NMF_H>0.5,canonical=="missense variant"))
syn_not5FU <- nrow(DATA_synonymous %>% filter(sampleid%in% nonFU, clonality=="SUBCLONAL",NMF_H>0.5))
nonsyn_not5FU <- nrow(DATA_nonsynonymous %>% filter(sampleid%in% nonFU, clonality=="SUBCLONAL",NMF_H>0.5,canonical=="missense variant"))

challenge.df = matrix(c(syn_5FU,syn_not5FU,nonsyn_5FU,nonsyn_not5FU), nrow = 2)
fisher.test(challenge.df, alternative ="greater")


library(dndscv)
DATA_synonymous %>% filter(sampleid%in% FU, NMF_breast_F>0.5) -> FU_induced_mutations_synonymous #NMF_colon_G
DATA_nonsynonymous %>% filter(sampleid%in% FU, NMF_breast_F>0.5) -> FU_induced_mutations_nonsynonymous #NMF_colon_G
df=rbind(FU_induced_mutations_nonsynonymous,FU_induced_mutations_synonymous)
head(df)
df <- df[,c("sampleid","chromosome","position", "ref_original","alt_original")]
colnames(df) <- c("sampleID", "chr", "pos", "ref", "mut")
df$chr <- substr(df$chr, 4, nchar(as.vector(df$chr)))
dndsout = dndscv(df)

sel_cv = dndsout$sel_cv
dndsout$globaldnds


sel_cv_5FU <- sel_cv
signif_genes_5FU = sel_cv_5FU[sel_cv_5FU$qglobal_cv<0.1, c("gene_name","qglobal_cv")]

