library(foreach)
library(MutationalPatterns)
library(BSgenome)
library(reshape)
library(gridExtra)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
bed_to_granges = function(bed_file)
{
  bed = read.table(bed_file, header = F, stringsAsFactors = F)
  chr = paste("chr", bed[,1], sep="")
  # Convert BED (0-based) start postion to Granges (1-based)
  start = bed[,2] + 1
  # In BED end position is excluded, in Granges end position is included -> +1 -1 -> no conversion needed
  end = bed[,3]
  GR = GRanges(chr, IRanges(start,end))  
  return(GR)
}

table_to_granges = function(df)
{
  table = df
  #chr = paste("chr", table[,1], sep="")
  # Convert BED (0-based) start postion to Granges (1-based)
  start = table[,2] + 1
  # In BED end position is excluded, in Granges end position is included -> +1 -1 -> no conversion needed
  end = table[,3]
  strand ="*"
  table[,4]
  
  #temp, GRanges(chromosome, IRanges(position, end = position), REF = ref, ALT 
  GR = GRanges(table[,1], IRanges(start,end),strand=strand,strand_info=as.factor(table[,4]))  
  return(GR)
}
strand_bias_test_pvalue = function(strand_occurrences)
{
  # These variables will be available at run-time, but not at compile-time.
  # To avoid compiling trouble, we initialize them to NULL.
  group = NULL
  type = NULL
  strand = NULL
  variable = NULL
  
  # statistical test for strand ratio
  # poisson test
  df_strand = reshape2::dcast(melt(strand_occurrences),
                              group + type ~ strand,
                              sum,
                              subset = plyr::.(variable == "no_mutations"))
  
  df_strand$total = df_strand[,3] + df_strand[,4]
  df_strand$ratio = df_strand[,3] / df_strand[,4]
  df_strand$p_poisson = apply(df_strand, 1, function(x) poisson.test(c(as.numeric(x[3]), as.numeric(x[4])), r=1)$p.value)
  df_strand$significant[df_strand$p_poisson < 0.01] = "*"
  df_strand$significant[df_strand$p_poisson >= 0.01] = " "
  
  return(df_strand)
}

#######
######
#extract CN>NT mutations
######

library(SPARQL, quietly = T)
library(stringr, quietly = T)
library(Biostrings, quietly = T)
library(dplyr,quietly = T)
library(foreach, quietly = T)

options(stringsAsFactors = F)

load("~/surfdrive/Shared/Sig17/HMF_data/cohort_analyse/somatic_clinical.RData")


endpoint <- "http://localhost:8890/sparql-auth"
auth_options <- curlOptions(userpwd="avanhoeck:XXXXXXXX")


#========= Query =========#



generate_TG_table <- function(cohort = Breast_patients){
  df_out=NULL
  df_out=data.frame()
  
  print(cohort)
  for (sampleID in unique(cohort)) {
    query <- sprintf(
      "
      PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
      PREFIX sampleid_pf: <http://sparqling-genomics/Sample/>
      PREFIX chromosome_pf: <http://rdf.biosemantics.org/data/genomeassemblies/hg19>
      PREFIX clonality_pf: <http://sparqling-genomics/Clonality/>
      PREFIX filter_pf: <http://sparqling-genomics/FilterItem/>
      PREFIX type_pf: <http://sparqling-genomics/VariantCallType/>
      PREFIX trinucleotidecontext_pf: <http://sparqling-genomics/NucleotideContext/>
      PREFIX seq_pf: <http://sparqling-genomics/Sequence/>
      
      SELECT 
      STRAFTER(STR(?sampleid), STR(sampleid_pf:)) AS ?sampleid
      STRAFTER(STR(?chromosome), \"http://rdf.biosemantics.org/data/genomeassemblies/hg19#\") AS ?chromosome
      ?position
      ?ref
      ?alt
      STRAFTER(STR(?clonality), \"http://sparqling-genomics/Clonality/\") AS ?clonality
      STRAFTER(STR(?trinucleotidecontext), \"http://sparqling-genomics/NucleotideContext/\") AS ?trinucleotidecontext
      STRAFTER(STR(?type), \"http://sparqling-genomics/VariantCallType/\") AS ?type
      STRAFTER(STR(?filter), \"http://sparqling-genomics/FilterItem/\") AS ?filter
      ?adjustedvaf
      
      FROM <http://hmfpatients/somaticvariant>
      WHERE 
      { ?row col:sampleid sampleid_pf:%s .
      ?row col:chromosome ?chromosome .
      ?row col:position ?position .
      ?row col:ref ?ref .
      ?row col:alt ?alt .
      ?row col:type ?type .
      ?row col:clonality ?clonality .
      ?row col:adjustedvaf ?adjustedvaf .
      ?row col:trinucleotidecontext ?trinucleotidecontext .
      ?row col:filter ?filter .
      BIND(sampleid_pf:%s AS ?sampleid)
      
      FILTER (regex(?trinucleotidecontext, \"C[A-Z]T\") || regex(?trinucleotidecontext, \"A[A-Z]G\" ))
      FILTER (?filter = filter_pf:PASS)
      FILTER (?type = type_pf:SNP)
      
      }
      ",sampleID,sampleID)
    query_out <- SPARQL(url = endpoint, curl_args = auth_options, query = query)$results
    df_out=rbind(df_out,query_out)
    strings <- c("T","C")
    reverse_strings <- c("A","G")
    df_out_1 <- df_out %>% filter(grepl('^C.T$', trinucleotidecontext) &  ref %in% strings)
    df_out_2 <- df_out %>% filter(grepl('^A.G$', trinucleotidecontext) &  ref %in% reverse_strings)
    
    df_out_final=rbind(df_out_1,df_out_2)
    
    
    
  }
  return(df_out_final)
  rm(df_out_1)
  rm(df_out_2)
}

Colon_T2G_all <- generate_TG_table(cohort = Colon_patients)
Breast_T2G_all <- generate_TG_table(cohort = Breast_patients)
Esophagus_T2G_all <- generate_TG_table(cohort = Esophagus_patients)

Colon_T2G_all <- Colon_T2G_all %>% filter(!grepl(',', alt) ) %>% as.data.frame()
Breast_T2G_all <- Breast_T2G_all %>% filter(!grepl(',', alt) ) %>% as.data.frame()
Esophagus_T2G_all <- Esophagus_T2G_all %>% filter(!grepl(',', alt) ) %>% as.data.frame()

vcf_list <- list.files("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/organoids/", pattern = ".vcf", full.names = TRUE)
vcf_files_names <- c("SI_5-FU_1","SI_5-FU_2")
vcfs_organoids <- read_vcfs_as_granges(vcf_list, vcf_files_names, genome = ref_genome)



#######
######
#extract CT>GT mutations
######
clonal <- data.frame(path=list.files("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/TG_CLONAL/", recursive = T, full.names = T))
subclonal <- data.frame(path=list.files("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/TG_SUBCLONAL/", recursive = T, full.names = T))
all_data=rbind(clonal,subclonal)
all_data$sampleId=gsub(pattern = "(.*R_)(.*)(_PASS_SNV.*)",replacement = "\\2",x = all_data$path)




genome_length <- sum(as.numeric(as.vector(as.data.frame(seqlengths(Hsapiens))[1:24,])))/1000000
somatic_clinical_colon$hypermutation = somatic_clinical_colon$mut_load/genome_length
somatic_clinical_subset = somatic_clinical_colon %>% 
  filter(Fluorouracil != "-2",hypermutation <=10)


Colon_patients = somatic_clinical_subset[which(somatic_clinical_subset$primaryTumorLocation == "Colon/Rectum" & somatic_clinical_subset$Signature.17 > 2000 & somatic_clinical_subset$Signature.17_rel > 0.25 & somatic_clinical_subset$Fluorouracil =="1" ),]$sampleId
Breast_patients = somatic_clinical_subset[which(somatic_clinical_subset$primaryTumorLocation == "Breast" & somatic_clinical_subset$Signature.17 > 2000 & somatic_clinical_subset$Signature.17_rel > 0.25 & somatic_clinical_subset$Fluorouracil =="1" ),]$sampleId
Esophagus_patients = somatic_clinical_subset[which(somatic_clinical_subset$primaryTumorLocation == "Esophagus" & somatic_clinical_subset$Signature.17 > 2000 & somatic_clinical_subset$Signature.17_rel > 0.25 & somatic_clinical_subset$Fluorouracil =="0" ),]$sampleId

Colon_all_data=as.character(all_data[which(all_data$sampleId %in% Colon_patients),]$path)
Breast_all_data=as.character(all_data[which(all_data$sampleId %in% Breast_patients),]$path)
Esophagus_all_data=as.character(all_data[which(all_data$sampleId %in% Esophagus_patients),]$path)

Colon_T2G <- foreach(f = c(Colon_all_data), .combine = 'rbind') %do% {
  df <- read.csv(f, sep= " ",colClasses=c(rep('character', 9)))
  return(df)
}
Breast_T2G <- foreach(f = c(Breast_all_data), .combine = 'rbind') %do% {
  df <- read.csv(f, sep= " ",colClasses=c(rep('character', 9)))
  return(df)
}
Esophagus_T2G <- foreach(f = c(Esophagus_all_data), .combine = 'rbind') %do% {
  df <- read.csv(f, sep= " ",colClasses=c(rep('character', 9)))
  return(df)
}

Colon_T2G=Colon_T2G[complete.cases(Colon_T2G), ]
Colon_T2G$position = as.numeric(Colon_T2G$position)
Breast_T2G=Breast_T2G[complete.cases(Breast_T2G), ]
Breast_T2G$position = as.numeric(Breast_T2G$position)
Esophagus_T2G=Esophagus_T2G[complete.cases(Esophagus_T2G), ]
Esophagus_T2G$position = as.numeric(Esophagus_T2G$position)


Colon_T2G_all 
Breast_T2G_all 
Esophagus_T2G_all 

#####
#####
#make GR objects per sampleId

Breast_T2G_IDnames = unique(Colon_T2G_all$sampleid)

GR <- GRangesList()
for (i in 1:length(Breast_T2G_IDnames)){
  temp <- Colon_T2G_all[Colon_T2G_all$sampleid==Breast_T2G_IDnames[i],]
  #print(head(temp))
  
  temp_GR = with(temp, GRanges(chromosome, IRanges(position, end = position), REF = ref, ALT = alt,clonality = clonality,trinucleotidecontext = trinucleotidecontext))
  seqlevelsStyle(temp_GR) <- "UCSC"
  auto <- extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
  temp_GR_all <- keepSeqlevels(temp_GR,auto, pruning.mode="coarse")
  #temp_GR_subclonal = temp_GR_all[temp_GR_all$clonality=="SUBCLONAL"]
  #temp_GR_clonal = temp_GR_all[temp_GR_all$clonality=="CLONAL"]
  
  temp_GR_all_GR_list = GRangesList(temp_GR_all)
  #name_granelist_all = paste(IDnames[i],"_all",sep = "")
  #names(temp_GR_all_GR_list) = name_granelist_all
  names(temp_GR_all_GR_list) = paste(Breast_T2G_IDnames[i],"_colon",sep = "")
  
  #temp_GR_subclonal_GR_list = GRangesList(temp_GR_subclonal)
  #name_granelist_subclonal = paste(IDnames[i],"_SUBCLONAL",sep = "")
  #names(temp_GR_subclonal_GR_list) = name_granelist_subclonal
  
  #temp_GR_clonal_GR_list = GRangesList(temp_GR_clonal)
  #name_granelist_clonal = paste(IDnames[i],"_CLONAL",sep = "")
  #names(temp_GR_clonal_GR_list) = name_granelist_clonal
  
  #GR <- c(GR,temp_GR_all_GR_list,temp_GR_subclonal_GR_list,temp_GR_clonal_GR_list)
  GR <- c(GR,temp_GR_all_GR_list)
}


GR_colon=GR
GR_breast=GR
GR_esophagus=GR
length(GR_colon)
length(GR_breast)
length(GR_esophagus)
GR <- c(GR_colon,GR_breast,GR_esophagus)


save(GR, file = "~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/genomic_range.RData")

#gragne objects
load("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/genomic_range.RData")


GR=GR_all
names(GR_all)
cancertype <- c(rep("colon", length(grep("colon", names(GR)))),
                rep("breast", length(grep("breast", names(GR)))),
                rep("esophagus", length(grep("esophagus", names(GR)))))

#Enrichment or depletion of mutations in genomic regions
###################################################
### code chunk number 94: download_using_biomaRt
###################################################
genome(GR) <- ref_genome

#test_bed_files: dit stuk script hieronder kan gebruikt worden om het script te normaliseren omdat het laden van de bed files heel lang duurt
surveyed_file <- bed_to_granges("~/Documents/ref_fasta/HG001_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel.bed")
seqlevelsStyle(surveyed_file) <- "UCSC"
genome(surveyed_file) <- ref_genome
surveyed_list <- rep(list(surveyed_file), length(GR))
surveyed_list = GRangesList(surveyed_list)


library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
GR_genes <- reduce(genes(txdb))
GR_cds <- reduce(cds(txdb))
GR_noncds <- setdiff(GR_genes,GR_cds)

gene_regions=GRangesList(GR_genes, GR_cds, GR_noncds)
names(gene_regions)=c("Genes","CDS", "non-CDS")

genome(gene_regions) <- ref_genome
genome(histo.regions) <- ref_genome

distr_gene_regions=genomic_distribution(GR,surveyed_list = surveyed_list, gene_regions)
distr_test_gene_regions=enrichment_depletion_test(distr_gene_regions, by=cancertype)
plot_enrichment_depletion(distr_test_gene_regions)

#other genomic context regions
#promotor and promotor flanking regions
promoter_g=readRDS(system.file("states/promoter_g_data.rds", package = "MutationalPatterns"))
flanking_g=readRDS(system.file("states/promoter_flanking_g_data.rds", package = "MutationalPatterns"))
promotor_regions=GRangesList(promoter_g, flanking_g)
names(promotor_regions)=c("Promoter", "Promoter flanking")
seqlevelsStyle(promotor_regions)="UCSC"
genome(promotor_regions) <- ref_genome
distr_promotor_regions=genomic_distribution(GR,surveyed_list = surveyed_list, promotor_regions)
distr_test_promotor_regions=enrichment_depletion_test(distr_promotor_regions, by=cancertype)
plot_enrichment_depletion(distr_test_promotor_regions)

#LAD regions
LAD = read.table("~/surfdrive/Shared/IPS_vs_Organoids/LAD/GSE22428_hESC.bed", header = FALSE,stringsAsFactors=FALSE) 
LAD.gr=makeGRangesFromDataFrame(LAD,keep.extra.columns=TRUE,
                                ignore.strand=TRUE, seqnames.field = 'V1',
                                start.field="V2",
                                end.field='V3')
LAD.regions = GRangesList(LAD.gr)
names(LAD.regions)=c("hESCs LADs")
plot_enrichment_depletion(distr_lad)
genome(LAD.regions) <- ref_genome
distr_LAD.regions=genomic_distribution(GR,surveyed_list = surveyed_list, LAD.regions)
distr_test_LAD.regions=enrichment_depletion_test(distr_LAD.regions, by=cancertype)
plot_enrichment_depletion(distr_test_LAD.regions)


#histones annotation (from Encode)
H3K4me1 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K4me1/H3K4me1.merged.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
H3K4me1.gr=makeGRangesFromDataFrame(H3K4me1,keep.extra.columns=FALSE,
                                    ignore.strand=TRUE, seqnames.field = 'V1',
                                    start.field="V2",
                                    end.field='V3')

H3K4me2 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K4me2/H3K4me2.merged.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
H3K4me2.gr=makeGRangesFromDataFrame(H3K4me2,keep.extra.columns=FALSE,
                                    ignore.strand=TRUE, seqnames.field = 'V1',
                                    start.field="V2",
                                    end.field='V3')

H3K4me3 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K4me3/H3K4me3.merged.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
H3K4me3.gr=makeGRangesFromDataFrame(H3K4me3,keep.extra.columns=FALSE,
                                    ignore.strand=TRUE, seqnames.field = 'V1',
                                    start.field="V2",
                                    end.field='V3')

H3K9me2 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K9me2/H3K9me2.merged.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
H3K9me2.gr=makeGRangesFromDataFrame(H3K9me2,keep.extra.columns=FALSE,
                                    ignore.strand=TRUE, seqnames.field = 'V1',
                                    start.field="V2",
                                    end.field='V3')

H3K27ac = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K27ac/H3K27ac.merged.bed", header = FALSE,
                     stringsAsFactors=FALSE) 
H3K27ac.gr=makeGRangesFromDataFrame(H3K27ac,keep.extra.columns=FALSE,
                                    ignore.strand=TRUE, seqnames.field = 'V1',
                                    start.field="V2",
                                    end.field='V3')

H3K27me3 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K27me3/H3K27me3.merged.bed", header = FALSE,
                      stringsAsFactors=FALSE) 
H3K27me3.gr=makeGRangesFromDataFrame(H3K27me3,keep.extra.columns=FALSE,
                                     ignore.strand=TRUE, seqnames.field = 'V1',
                                     start.field="V2",
                                     end.field='V3')

H3K36me3 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K36me3/H3K36me3.merged.bed", header = FALSE,
                      stringsAsFactors=FALSE) 
H3K36me3.gr=makeGRangesFromDataFrame(H3K36me3,keep.extra.columns=FALSE,
                                     ignore.strand=TRUE, seqnames.field = 'V1',
                                     start.field="V2",
                                     end.field='V3')

H3K79me2 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3K79me2/H3K79me2.merged.bed", header = FALSE,
                      stringsAsFactors=FALSE) 
H3K79me2.gr=makeGRangesFromDataFrame(H3K79me2,keep.extra.columns=FALSE,
                                     ignore.strand=TRUE, seqnames.field = 'V1',
                                     start.field="V2",
                                     end.field='V3')

H4K20me1 = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H4K20me1/H4K20me1.merged.bed", header = FALSE,
                      stringsAsFactors=FALSE) 
H4K20me1.gr=makeGRangesFromDataFrame(H4K20me1,keep.extra.columns=FALSE,
                                     ignore.strand=TRUE, seqnames.field = 'V1',
                                     start.field="V2",
                                     end.field='V3')

H2AFZ = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H2AFZ/H2AFZ.merged.bed", header = FALSE,
                   stringsAsFactors=FALSE) 
H2AFZ.gr=makeGRangesFromDataFrame(H2AFZ,keep.extra.columns=FALSE,
                                  ignore.strand=TRUE, seqnames.field = 'V1',
                                  start.field="V2",
                                  end.field='V3')

H3F3A = read.table("~/surfdrive/Shared/IPS_vs_Organoids/encode/H3F3A/H3F3A.merged.bed", header = FALSE,
                   stringsAsFactors=FALSE) 
H3F3A.gr=makeGRangesFromDataFrame(H3F3A,keep.extra.columns=FALSE,
                                  ignore.strand=TRUE, seqnames.field = 'V1',
                                  start.field="V2",
                                  end.field='V3')


histo.regions = GRangesList(H3K4me1.gr, H3K4me2.gr, H3K4me3.gr, H3K9me2.gr, H3K27ac.gr, H3K27me3.gr, H3K36me3.gr, H3K79me2.gr, H4K20me1.gr, H2AFZ.gr, H3F3A.gr)
names(histo.regions)=c("H3K4me1","H3K4me2", "H3K4me3", "H3K9me2", "H3K27ac", "H3K27me3", "H3K36me3","H3K79me2", "H4K20me1", "H2AFZ", "H3F3A")
genome(histo.regions) <- ref_genome
#2A
#TFs annotation (from Encode) + DHS from chip-atlas
genome(histo.regions) <- ref_genome
distr_histo.regions=genomic_distribution(GR,surveyed_list = surveyed_list, histo.regions)
distr_test_histo.regions=enrichment_depletion_test(distr_histo.regions, by=cancertype)
plot_enrichment_depletion(distr_test_histo.regions)


RepliSeq_early = bed_to_granges("~/Documents/ref_fasta/all_RepliSeq_median_early_merged.bed")
RepliSeq_intermediate = bed_to_granges("~/Documents/ref_fasta/all_RepliSeq_median_intermediate_merged.bed")
RepliSeq_late = bed_to_granges("~/Documents/ref_fasta/all_RepliSeq_median_late_merged.bed")

RepliSeq = GRangesList(early=RepliSeq_early,intermediate=RepliSeq_intermediate,late=RepliSeq_late)
auto <- extractSeqlevelsByGroup(species="Homo_sapiens",style="UCSC",group="auto")
RepliSeq <- lapply(RepliSeq, function(x) keepSeqlevels(x, auto, pruning.mode="coarse"))
RepliSeq_grl=GRangesList(RepliSeq)
genome(RepliSeq_grl)=ref_genome
distr_RepliSeq_grl=genomic_distribution(GR,surveyed_list = surveyed_list, RepliSeq_grl)
distr_test_RepliSeq_grl=enrichment_depletion_test(distr_RepliSeq_grl, by=cancertype)
plot_enrichment_depletion(distr_test_RepliSeq_grl)


# ------- REPLICATION STRAND BIAS -------

###replicseq data
repli_file = system.file("extdata/ReplicationDirectionRegions.bed",package = "MutationalPatterns")
repli_strand = read.table(repli_file, header = TRUE)
# Store in GRanges object
repli_strand_granges = GRanges(seqnames = repli_strand$Chr,
                               ranges = IRanges(start = repli_strand$Start + 1,end = repli_strand$Stop),
                               strand_info = as.factor(repli_strand$Class))



###tonkova et al
Haradhvala <- read.table("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/bsblabludwig-replicationasymmetry-e89d5d83f084/data/tableTerritories_Besnard1k_territories_1000_bins.txt", header = TRUE)

Haradhvala %>% dplyr::filter(abs(nBinsFromORI) <= 5 ) %>% 
  dplyr::select(tChr,tPos0,tPos1,tIsLeft,tIsRight,nBinsFromORI) %>% 
  dplyr::mutate(sum_left_rigt = tIsLeft + tIsRight)%>% 
  dplyr::filter(sum_left_rigt==1 ) %>%
  dplyr::mutate(strand_info = ifelse(tIsLeft == "1", "left", ifelse(tIsRight == "1", "right", "NA"))) %>% 
  dplyr::select(tChr,tPos0,tPos1,strand_info) %>% 
  as.data.frame() -> 
  Haradhvala_repli_strand



Haradhvala_repli_strand_granges=table_to_granges(Haradhvala_repli_strand)

#mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))

# UCSC seqlevelsstyle
repli_strand_granges=Haradhvala_repli_strand_granges
seqlevelsStyle(repli_strand_granges) = "UCSC"

# Make mutation count matrix with transcriptional strand information 
# (96 trinucleotides * 2 strands = 192 features).
mut_mat_s_rep <- mut_matrix_stranded(GR, ref_genome, repli_strand_granges,mode = "replication")
#mut_mat_s_rep_all=mut_mat_s_rep
#mut_mat_s_rep_all[,84]
#mut_mat_s_rep[,84]
#mut_mat_s_rep[c("C[T>G]T-left"), 1:8]
#mut_mat_s_rep[c("C[T>G]T-right"), 1:8]

# Count the number of mutations on each strand, per tissue, per mutation type
strand_counts_rep <- strand_occurrences(mut_mat_s_rep, by=cancertype)
strand_counts_rep$type <- sprintf("C[%s]T", strand_counts_rep$type)
# Perform Poisson test for strand asymmetry significance testing:
#strand_bias_rep <- strand_bias_test_pvalue(strand_counts_rep)
strand_bias_rep <- strand_bias_test(strand_counts_rep)

# Plot the mutation spectrum with strand distinction:
ps1 <- plot_strand(strand_counts_rep, mode = "relative")
#Plot the effect size (log2(untranscribed/transcribed) of the strand bias. Asteriks indicate significant strand bias.
ps2 <- plot_strand_bias(strand_bias_rep)
# Combine the plots into one figure:
grid.arrange(ps1, ps2)

pdf(paste("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/plots/", "rep_strand_biasTNtoNGmuts.pdf", sep = ""), width = 8, height = 7, useDingbats = F)
grid.arrange(ps1, ps2)
dev.off()



# ------- Transcription STRAND BIAS -------

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
genome(genes_hg19) <- ref_genome
genome(vcfs_all) <- ref_genome
mut_mat_s <- mut_matrix_stranded(GR, ref_genome, genes_hg19)
mut_mat_s
strand_counts <- strand_occurrences(mut_mat_s, by=cancertype)
strand_counts$type <- sprintf("C[%s]T", strand_counts$type)

strand_bias <- strand_bias_test_pvalue(strand_counts)
strand_bias
ps1 <- plot_strand(strand_counts, mode = "relative")
ps2 <- plot_strand_bias(strand_bias)
grid.arrange(ps1, ps2)
pdf(paste("~/surfdrive/Shared/Sig17/HMF_data/Compare_signatures/plots/", "trans_strand_bias_TNtoNGmuts.pdf", sep = ""), width = 8, height = 7, useDingbats = F)
grid.arrange(ps1, ps2)
dev.off()

