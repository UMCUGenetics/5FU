#!/usr/bin/env Rscript
library(GenomicRanges)
library(VariantAnnotation)
library(ggplot2)
library(reshape2)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)
library("readxl")
library(dplyr)


###output_dir
dirpath="~/surfdrive/Shared/Sig17/HMF_data/somatic_organoid_vcfs/HMF_pipeline/plots/March/"

#load COSMIC signatures
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

#load de novo obtained signatures
cancer_signatures_denovo=as.matrix(as.data.frame(read.table(file =  "~/surfdrive/Shared/Sig17/HMF_data/extract_signatures/colon_breast_invitro/nmf_rank3_16_signature.txt",sep = "\t", header = TRUE, row.names = 1)))

###Functions
cos_sim = function(x, y){
  res = x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res = as.numeric(res)
  return(res)
}
colorpalette = c("Signature.1" =  '#8dd3c7',
                 "Signature.2" =  '#ffffb3',
                 "Signature.3" =  '#bebada',
                 "Signature.4" =  '#fb8072',
                 "Signature.5" =  '#80b1d3',
                 "Signature.6" =  '#fdb462',
                 "Signature.7" =  '#b3de69',
                 "Signature.8" =  '#fccde5',
                 "Signature.9" =  '#d9d9d9',
                 "Signature.10" = '#ff1417' ,
                 "Signature.11" = '#ff6611' ,
                 "Signature.12" = '#c4ff00' ,
                 "Signature.13" = '#ff8844' ,
                 "Signature.14" = '#ffee55' ,
                 "Signature.15" = '#ffff99' ,
                 "Signature.16" = '#78FA37' ,
                 "Signature.17" = '#aacc22' ,
                 "Signature.18" = '#bbdd77' ,
                 "Signature.19" = '#c8cf82' ,
                 "Signature.20" = '#92a77e' ,
                 "Signature.21" = '#5599ee' ,
                 "Signature.22" = '#0088cc' ,
                 "Signature.23" = '#226688' ,
                 "Signature.24" = '#175279' ,
                 "Signature.25" = '#557777' ,
                 "Signature.26" = '#ddbb33' ,
                 "Signature.27" = '#d3a76d' ,
                 "Signature.28" = '#a9834b' ,
                 "Signature.29" = '#aa6688',
                 "Signature.30" = '#767676',
                 "Signature.A" = '#458B00' ,
                 "Signature.B" = '#D2691E' ,
                 "Signature.C" = '#6495ED' ,
                 "Signature.D" = '#A2CD5A' ,
                 "Signature.E" = '#CD3333' ,
                 "Signature.F" = '#7AC5CD' ,
                 "Signature.G" = '#009ACD' ,
                 "Signature.H" = '#CD2626' ,
                 "Signature.I" = '#FFB90F' ,
                 "Signature.J" = '#76EEC6' ,
                 "Signature.K" = '#EEB422' ,
                 "Signature.L" = '#97FFFF' ,
                 "Signature.M" = '#E9967A' ,
                 "Signature.N" = '#5F9EA0')


#import SNVs
vcf_organoids <- list.files("~/surfdrive/Shared/Sig17/HMF_data/somatic_organoid_vcfs/filtered_pass/invitro/organoid_VCFs/", pattern = ".vcf", full.names = TRUE)
vcf_files_names <- substr(basename(vcf_list), 1, nchar(basename(vcf_list)) - 4) 
vcf_files_names <- sub("_pass_filtered_SNV_VAF30_70","",vcf_files_names)
vcfs_SC <- read_vcfs_as_granges(vcf_list, vcf_files_names, genome = ref_genome)
treatment <- c(rep("Control", 6),rep("5-FU", 2))

#plot type occurences
type_occurrences <- mut_type_occurrences(vcfs_SC, ref_genome)
type_occurrences
plot_spectrum(type_occurrences)
plot_spectrum(type_occurrences, CT = TRUE, legend = FALSE)

#total number of SNVs
length(vcfs_SC$"STE072-control-p17_5-FU-2-625-8")
length(vcfs_SC$"STE072-control-p17_5-FU-3-625-7")

#create 96-matrix
auto <- extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
vcfs <- lapply(vcfs_SC, function(x) keepSeqlevels(x, auto, pruning.mode="coarse"))
vcfs_mm <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
colnames(vcfs_mm)
colSums(vcfs_mm)
plot_96_profile(vcfs_mm)
cos_sim(vcfs_mm[,7], vcfs_mm[,8])

vcfs_mm_df <- as.data.frame(vcfs_mm)
vcfs_mm_df <- cbind(as.data.frame(rowSums(vcfs_mm_df[grep("5-FU", names(vcfs_mm_df),value = T)])),
                    as.data.frame(rowSums(vcfs_mm_df[grep("^STE00", names(vcfs_mm_df),value = T)])))
colnames(vcfs_mm_df) <- c("5-FU", "Control")
plot_96_profile(vcfs_mm_df)
cos_sim(vcfs_mm_df[,1], vcfs_mm_df[,2])
plot_compare_profiles(vcfs_mm_df[,1], 
                      vcfs_mm_df[,2], 
                      profile_names = c("5-FU", "Control"),
                      condensed = TRUE)

signature_5FU <- vcfs_mm_df[,1]
signature_control <- vcfs_mm_df[,2]
relative_5FU = signature_5FU / sum(signature_5FU)
relative_control = signature_control / sum(signature_control)
diff = relative_5FU - relative_control
diff_neg <- diff
diff_neg[diff_neg>0] <- 0
diff_neg <- abs(diff_neg)
diff_pos <- diff
diff_pos[diff_pos<0] <- 0
format(diff_pos, scientific = F, digits = 3)
vcfs_mm_df$diff=diff
vcfs_mm_df$diff_5FU=diff_pos
vcfs_mm_df$diff_control=diff_neg

#plot in vito, in vitro and cosmic 17 signatures
signatures_bucket=data.frame(InVitro_5FU=diff_pos,
                             InVivo_5FU=as.data.frame(cancer_signatures_denovo)$NMF_H,
                             COSMIC_17=as.data.frame(cancer_signatures)$Signature.17)
plot_96_profile(as.matrix(signatures_bucket), condensed = T)


#compare in-vitro 5FU signature to COSMIC and de novo obtained signatures
cos_sim(signatures_bucket$InVitro_5FU, signatures_bucket$COSMIC_17)
cos_sim(signatures_bucket$InVitro_5FU, signatures_bucket$InVivo_5FU)
#pearson
cor(signatures_bucket$InVitro_5FU, signatures_bucket$COSMIC_17, method = c("pearson"))
cor(signatures_bucket$InVitro_5FU, signatures_bucket$InVivo_5FU, method=c("pearson"))



norm_mut_matrix_new <- as.data.frame(vcfs_mm)
norm_mut_matrix_new <- cbind(as.data.frame(rowSums(norm_mut_matrix_new[grep("5-FU", names(norm_mut_matrix_new),value = T)])),
                    as.data.frame(rowSums(norm_mut_matrix_new[grep("^STE00", names(norm_mut_matrix_new),value = T)])))
colnames(norm_mut_matrix_new) <- c("5-FU", "Control")
norm_mut_matrix_new$diff_pos = diff_pos
norm_mut_matrix_new$diff_neg = diff_neg
norm_mut_matrix_new = as.data.frame(apply(as.matrix(norm_mut_matrix_new), 2, function(x) x / sum(x) ))
norm_mut_matrix_new$diff_neg=norm_mut_matrix_new$diff_neg*(-1)
norm_mut_matrix_new$diff=norm_mut_matrix_new$diff_pos+norm_mut_matrix_new$diff_neg
norm_mut_matrix_new$diff_pos=NULL
norm_mut_matrix_new$diff_neg=NULL


colnames(norm_mut_matrix_new) = c("5-FU", "Control", "Difference")
norm_mut_matrix_new
colors=COLORS6
context = CONTEXTS_96
substitution = rep(SUBSTITUTIONS, each=16)
substring(context, 2, 2) = "."
df = data.frame(substitution = substitution, context = context)
rownames(norm_mut_matrix_new) = NULL
df2 = cbind(df, as.data.frame(norm_mut_matrix_new))
df3 = melt(df2, id.vars = c("substitution", "context"))

ymax=0.4
y.upper.limit <- diff(range(df3$value)) * 0.05 + max(df3$value)
y.lower.limit <- 0 - diff(range(df3$value)) * 0.2


plot_signatures = ggplot(data=df3, aes(x=context,
                                       y=value,
                                       fill=substitution,
                                       width=1)) +
  geom_bar(stat="identity", colour="black", size=.2) +
  #geom_col(position = "dodge", colour = "black") +
  scale_fill_manual(values=colors) +
  facet_grid(variable ~ substitution) +
  ylab("Relative contribution") +
  #coord_cartesian(ylim=c(0,ymax)) +
  coord_cartesian(ylim = c(y.lower.limit, y.upper.limit)) +
  scale_y_continuous(breaks=seq(0, ymax,0.1),labels = scales::percent_format(accuracy = 1)) +
  #scale_y_continuous(breaks = pretty(dat$y, n = 10)) +
  #geom_hline(yintercept = y.upper.limit) +
  
  # no legend
  guides(fill=FALSE) +
  # white background
  theme_bw() +
  # format text
  theme(axis.title.y=element_text(size=70,vjust=3),
        axis.text.y=element_text(size=70),
        axis.title.x=element_text(size=70),
        axis.text.x = element_text(colour="grey20",size=25,angle=45,hjust=.5,vjust=.5,face="plain"),
        axis.ticks.y = element_line(size=2),
        strip.text.x=element_text(size=50),
        strip.text.y=element_text(size=50,vjust=2),
        panel.grid.major.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(4, "lines"),
        axis.line = element_line(color = 'grey20',size = 2),
        strip.background=element_blank())

plot_name="invitro_signatures_figure1"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 50, height = 30) 
plot(plot_signatures)
dev.off()


cos_sim_samples_signatures = cos_sim_matrix(cancer_signatures_denovo,as.matrix(signatures_bucket$InVitro_5FU))
colnames(cos_sim_samples_signatures) <- c("  5-FU in vitro  ")
cos_sim_samples_signatures.m <- melt(cos_sim_samples_signatures)
cos_sim_samples_signatures.m
cos_sim_samples_signatures.m$Var2 = factor(cos_sim_samples_signatures.m$Var2, levels = c("  5-FU in vitro  "))
heatmap = ggplot(cos_sim_samples_signatures.m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1.01),na.value = "white") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x=NULL, y=NULL)+
  geom_text(aes(label = round(value, 2)), size = 4.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1.05,vjust = 1.1,size=12),
        axis.text.y = element_text(angle = 0, hjust = 0,size=20),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_name="cosine_sim_Invitro_vs_denovo"
pdf(sprintf("%s%s_2.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 1.2) 
plot(heatmap)
dev.off()


###########
##########

#Cosmic
cos_sim(signatures_bucket$InVivo_5FU, signatures_bucket$COSMIC_17)
#pearson
cor(signatures_bucket$InVivo_5FU, signatures_bucket$COSMIC_17, method = c("pearson"))


###################
###################
#compare signatures

#compare in vitro and in vivo signature to COSMICs
cos_sim_samples_signatures = cos_sim_matrix(cancer_signatures,as.matrix(signatures_bucket[,c("InVitro_5FU","InVivo_5FU")]))
colnames(cos_sim_samples_signatures) <- c("  5-FU in vitro  ","  5-FU in vivo  ")
cos_sim_samples_signatures.m <- melt(cos_sim_samples_signatures)
cos_sim_samples_signatures.m
cos_sim_samples_signatures.m$Var2 = factor(cos_sim_samples_signatures.m$Var2, levels = c("  5-FU in vivo  ","  5-FU in vitro  "))

heatmap = ggplot(cos_sim_samples_signatures.m, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "YlGnBu", direction = 1, name = "Cosine \nsimilarity", limits = c(0,1.01),na.value = "white") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x=NULL, y=NULL)+
  geom_text(aes(label = round(value, 2)), size = 4.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1.05,vjust = 1.1,size=12),
        axis.text.y = element_text(angle = 0, hjust = 0,size=20),
        panel.border = element_blank(),panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_name="cosine_sim_invitro_invivo_vs_COSMIC"
pdf(sprintf("%s%s.pdf",dirpath,plot_name) , useDingbats = F, width = 10, height = 1.2) 
plot(heatmap)
dev.off()
