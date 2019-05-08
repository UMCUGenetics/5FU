library(MutationalPatterns)
library(GenomicRanges)
library(ggplot2) 
library(gridExtra)
library(grid)
library(dplyr)
library("GenomeInfoDb")
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
#############
#denovosignature analysis
#############


dirpath="~/surfdrive/Shared/Sig17/HMF_data/extract_signatures/colon_breast_invitro2/"


########################
########################
####BREAST + colon + invitro
########################
########################
load("~/surfdrive/Shared/Sig17/HMF_data/cohort_analyse/somatic_clinical.RData")

#genome_length <- 2858674662/1000000
#somatic_clinical$hypermutation = somatic_clinical$mut_load/genome_length
somatic_clinical_subset = somatic_clinical %>% 
  dplyr::filter(primaryTumorLocation == "Breast"|primaryTumorLocation == "Colon/Rectum",Fluorouracil != "-2",Cohort=="Metastatic")
names(somatic_clinical)
somatic_clinical_subset %>% group_by(Fluorouracil) %>%
  summarize(count = n(), 
            mean_sign = mean(NMF_H),median_sign = median(NMF_H), SD_sign = sd(NMF_H),
            sum_sign=sum(NMF_H),max_sign=max(NMF_H),min_sign=min(NMF_H),
            mean_mut_load = mean(mut_load), SD_mut_load = sd(mut_load),
            median_mut_load = median(mut_load),
            sum_mut_load=sum(mut_load))

######################
######################


temp = somatic_clinical_subset

temp=cbind(temp[grep(">", names(temp),value = T)],temp[c("sampleId","primaryTumorLocation","Cohort")])
temp=cbind(temp[,-grep("_rel", names(temp))])
temp=temp[temp$Cohort=="Metastatic",]
temp$primaryTumorLocation

##denovo mutations
mut_mat <- as.matrix(temp[grep(">", names(temp),value = T)]) 
row.names(mut_mat) =temp$sampleId

#add invitro organoid data
vcf_organoids <- list.files("~/surfdrive/Shared/Sig17/HMF_data/somatic_organoid_vcfs/filtered_pass/invitro/organoid_VCFs/", pattern = ".vcf", full.names = TRUE)
vcf_files_names <- substr(basename(vcf_list), 1, nchar(basename(vcf_list)) - 4) 
vcf_files_names <- sub("_pass_filtered_SNV_VAF30_70","",vcf_files_names)
vcfs_SC <- read_vcfs_as_granges(vcf_list, vcf_files_names, genome = ref_genome)
treatment <- c(rep("Control", 6),rep("5-FU", 2))
auto <- extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
vcf_organoids <- lapply(vcf_organoids, function(x) keepSeqlevels(x, auto, pruning.mode="coarse"))
vcf_organoids_matrix <- t(mut_matrix(vcf_list = vcf_organoids, ref_genome = ref_genome))

mut_mat_invivo_invitro = rbind(mut_mat,vcf_organoids_matrix)

mut_mat <- as.matrix(read.table(file = "~/surfdrive/Shared/Sig17/HMF_data/matrices/colon_breast_invitro_mut_matrix",sep = "\t", header = TRUE, row.names = 1))
mut_mat=mut_mat_invivo_invitro+ 0.0001

library("NMF")
library("gridExtra")

rank_numbers=c("6","7","8","9","10","11","12","13","14","15","16","17","18","19")

for(i in 1:length(rank_numbers)){
  plotname="nmf_rank3_"
  rank_number=as.numeric(rank_numbers[i])
  print("start running nmf with rank number")
  print(rank_number)
  res <- nmf(t(mut_mat), rank=rank_number, method="brunet", nrun=100, seed=123456)
  print("nmf ended")
  signatures_res <- round(NMF::basis(res))
  contribution_res <- NMF::coef(res)
  reconstructed_res <- signatures_res %*% contribution_res
  nmf_res <- list(signatures = signatures_res,
                  contribution = contribution_res,
                  reconstructed = reconstructed_res)
  colnames(nmf_res$signatures) <- paste("NMF_", toupper(letters[1:rank_number]), sep="")
  signatures_res_rel = apply(signatures_res, 2, function(x) x / sum(x))
  colnames(signatures_res_rel) <- paste("NMF_", toupper(letters[1:rank_number]), sep="")
  contribution_res_abs <- round(coef(res)*colSums(basis(res)))
  reconstructed_res_abs <- signatures_res_rel %*% contribution_res_abs
  row.names(contribution_res_abs)=colnames(signatures_res_rel)
  nmf_res_abs <- list(signatures = signatures_res_rel,
                      contribution = contribution_res_abs,
                      reconstructed = reconstructed_res_abs)
  plot_96=plot_96_profile(nmf_res_abs$signatures,condensed = T)
  pdf(sprintf("%s%s%s_signature_contribution_plot.pdf",dirpath,plotname,rank_number) , useDingbats = F, width = 18, height = 13) 
  grid.draw(plot_96)
  dev.off()
  
  plot_contr=plot_contribution(nmf_res_abs$contribution, nmf_res_abs$signature,mode = "absolute")
  pdf(sprintf("%s%s%s_mutation_contribution_plot.pdf",dirpath,plotname,rank_number) , useDingbats = F, width = 20, height = 8) 
  grid.draw(plot_contr)
  dev.off()
  
  write.table(as.data.frame(nmf_res_abs$signatures), file = sprintf("%s%s%s_signature.txt",dirpath,plotname,rank_number),sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)
  write.table(as.data.frame(nmf_res_abs$contribution), file = sprintf("%s%s%s_contribution.txt",dirpath,plotname,rank_number),sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)

  cos_sim_all = round(cos_sim_matrix(nmf_res_abs$signatures, cancer_signatures),3)
  write.table(as.data.frame(cos_sim_all), file = sprintf("%s%s%s_cos_sim_cosmic.txt",dirpath,plotname,rank_number),sep = "\t", col.names = NA,qmethod = "double", quote = FALSE)
  cossim=plot_cosine_heatmap(cos_sim_all)
  pdf(sprintf("%s%s%s_cosine_sim_plot.pdf",dirpath,plotname,rank_number) , useDingbats = F, width = 8, height = 5) 
  grid.draw(cossim)
  dev.off()
}

  
