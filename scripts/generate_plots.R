library(reshape2)
library(plyr)
library(ggpubr)
library(dplyr)
library(BSgenome)
library(ggplot2)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

library(nlme)
library(RColorBrewer)



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
######################
######################
CheckTMB <- function(data = somatic_clinical_subset, pdf.path = NULL){
  
  data <- data %>% filter(Cohort == "Metastatic",hasSystemicPreTreatment=="Yes",TMB <=10) %>% dplyr::select(Fluorouracil,NMF_A:NMF_P) #NMF_colon_A:NMF_colon_J  #NMF_breast_A:NMF_breast_J

  data_log <- cbind(data[1],log10(data[-1]+1))
  print("T-test p values")
  ttestresults <- as.list(lapply(names(data)[-1],function(x)
    t.test(as.formula(paste(x,"Fluorouracil",sep="~")),data=data_log)))
  ttest.pval <- sapply(ttestresults, '[[', 'p.value')
  print(rbind(names(data)[-1],ttest.pval))
  
  print("Wilcox test p values")
  wilcoxresults <- as.list(lapply(names(data_log)[-1],function(x)
    wilcox.test(as.formula(paste(x,"Fluorouracil",sep="~")),data=data_log,alternative = "two.sided")))
  wilcox.test.pval <- sapply(wilcoxresults, '[[', 'p.value')
  print(rbind(names(data)[-1],wilcox.test.pval))

  
  
  checkTMB_m = melt(data,id=c("Fluorouracil"))
  checkTMB_m$sign_contribution=log10(checkTMB_m$value+1)
  checkTMB_m <- checkTMB_m %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  checkTMB_m$Fluorouracil_new <- factor(checkTMB_m$Fluorouracil_new, levels = c("Not 5-FU pretreated","5-FU pretreated"))
  
  plot <- ggboxplot(checkTMB_m, x = "Fluorouracil_new", y = "sign_contribution",
            color = "black",fill = "Fluorouracil_new",facet.by="variable")+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
    #stat_compare_means(aes(label = ..p.signif..))+
    stat_compare_means(method = "wilcox.test")+
    xlab(c(""))+
    ylab(c("Signature mutation contribution"))+
    #scale_y_continuous(limits = c(0, 0.5))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,8,10)
    plot(plot)
    dev.off()
  }
  
}
CheckTMB(data = somatic_clinical_subset)
CheckTMB(data = somatic_clinical_subset,pdf.path=sprintf("%s%s_CheckTMB.pdf",dirpath,plot_name)) 

plotting_data_cohort=data.frame(sampleId=somatic_clinical_subset$sampleId,
                                SIGN_5FU_denovo=somatic_clinical_subset$NMF_H,    #breast: NMF_breast_F    colon: NMF_colon_C #NMF_H
                                SIGN_5FU_denovo_rel=somatic_clinical_subset$NMF_H_rel,
                                SIGN_17=somatic_clinical_subset$Signature.17,
                                SIGN_17_rel=somatic_clinical_subset$Signature.17_rel,
                                T2G=somatic_clinical_subset$`C[T>G]T`,
                                T2G_rel=somatic_clinical_subset$`C[T>G]T_rel`,
                                Fluorouracil=as.factor(somatic_clinical_subset$Fluorouracil),
                                Fluorouracil_nr_predrugs=somatic_clinical_subset$Fluorouracil_nr_predrugs,
                                Fluorouracil_predrugs_exposure_days=somatic_clinical_subset$Fluorouracil_predrugs_exposure_days,
                                mut_load=somatic_clinical_subset$mut_load,
                                Capecitabine=somatic_clinical_subset$Capecitabine_subset,
                                Fluoruracil=somatic_clinical_subset$Fluorouracil_subset,
                                Folinic_acid=somatic_clinical_subset$FolinicUacid_subset,
                                copyNumber_MTHFR=somatic_clinical_subset$copyNumber_MTHFR,
                                copyNumber_DPYD=somatic_clinical_subset$copyNumber_DPYD,
                                copyNumber_TYMP=somatic_clinical_subset$copyNumber_TYMP,
                                copyNumber_TK1=somatic_clinical_subset$copyNumber_TK1,
                                copyNumber_DTYMK=somatic_clinical_subset$copyNumber_DTYMK,
                                copyNumber_TYMS=somatic_clinical_subset$copyNumber_TYMS,
                                Cohort=somatic_clinical_subset$Cohort,
                                TMB=somatic_clinical_subset$TMB,
                                hasSystemicPreTreatment=somatic_clinical_subset$hasSystemicPreTreatment)


#main figures:
plotMutLoadAnd5FuContrib <- function(data = plotting_data_cohort, pdf.path = NULL){
  
  data <- data %>% 
    filter(Cohort == "Metastatic") %>% arrange(desc(SIGN_5FU_denovo))
  data$sampleId <- factor(data$sampleId, levels = as.character(data$sampleId))
  
  data$log_mut_load <- log10(data$mut_load)
  data$SIGN_5FU_denovo_magn <- data$SIGN_5FU_denovo/1000
  
  scale_factor <- max(data$SIGN_5FU_denovo_magn) / max(data$log_mut_load)
  
  colors <- list(
    bar.5fu.treated = '#2E6DAE',
    bar.non.treated = '#E6F5C1',
    line.outline = '#49B7C3',
    line.fill = '#6DC4BE'
  )
  
  ## Translate vector of 1/0 to vector of colors
  bar_colors <- ifelse(
    as.integer(as.character(data$Fluorouracil)),
    colors$bar.5fu.treated,
    colors$bar.non.treated
  )
  
  plot <- ggplot(data, aes(x=sampleId)) + 
    
    ## Bars: mut load; 5FU treated yes/no
    geom_bar(
      aes(y=log_mut_load, fill=Fluorouracil), stat='identity', 
      width=1, position = position_nudge(x = 0.5), ## Align bars to the right of tick
      fill=bar_colors
    ) +
    ylab( expression(Bars:~log[10]~(mutational~load)) ) +
    #ylab('') +
    xlab('') +
    
    ## Line/ribbon overlay: mut load
    ## Scale 5FU contrib down while making primary axis, then back up while making secondary axis
    ## so that axes are scaled to maximum
    geom_ribbon(
      aes(ymin=0, ymax=SIGN_5FU_denovo_magn/scale_factor, group = 1),
      fill=colors$line.fill
    ) +
    geom_line(
      aes(y=SIGN_5FU_denovo_magn/scale_factor, group = 1),
      color=colors$line.outline, size=1
    ) +
    scale_y_continuous(
      position = 'right',
      sec.axis=sec_axis(
        ~.*scale_factor,
        name=expression(Line:~5-FU~signature~contribution~(X10^{3}))),
        #name=expression(Bars:~log[10]~(mutational~load))),
      expand=c(0,0),
      limits = c(0, 5)
    ) +
    
    theme(
      legend.position = 'none',
      axis.text.x=element_blank(),
      axis.text.y.left = element_text(size = 10),
      axis.text.y.right = element_text(size = 10),
      axis.ticks.x=element_blank(),
      panel.background=element_blank(),
      panel.border=element_rect(fill=NA),
      axis.title.y.right=element_text(angle=90)
    )
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,8,4)
    plot(plot)
    dev.off()
  }
  
}
plotMutLoadAnd5FuContrib(data = plotting_data_cohort) 
plotMutLoadAnd5FuContrib(data = plotting_data_cohort,pdf.path=sprintf("%s%s_plotMutLoadAnd5FuContrib.pdf",dirpath,plot_name)) 


plot5Fucomparison <- function(data = plotting_data_cohort,mode="absolute",control_cohort="No", pdf.path = NULL){
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  #check number of samples
  
  
  
  if (mode=="relative")
  {
    #select_test
    #Kruskal???Wallis
    #res.dunn <- dunnTest(SIGN_5FU_denovo_rel ~ Fluorouracil, data = data, method = "holm")
    #control_compare <- which(grepl("0 -",res.dunn$res$Comparison))
    #res <- data.frame(res.dunn$res[control_compare,])
    
    #wilcox
    #wilcox <- wilcox.test(SIGN_5FU_denovo_rel ~ Fluorouracil, data = plotting_data_cohort[which(plotting_data_cohort$Cohort=="Metastatic"),])
    #print(wilcox$p.value)
    
    
    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo_rel ~ Fluorouracil_new, data = data));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    nn = data %>% group_by(Fluorouracil_new) %>% tally()
    
    set_max <- data %>% group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_rel),median_sign = median(SIGN_5FU_denovo_rel), SD_sign = sd(SIGN_5FU_denovo_rel),max(SIGN_5FU_denovo_rel), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_rel)`)
    
    if(number<4){
      Y_pos = c(98/100,89/100,84/100)
      #Y_pos = set_max*Y_pos
      Y_pos = 1*Y_pos
      Y_pos = Y_pos[1:number]
    }
    

    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "SIGN_5FU_denovo_rel",
              color = "black",fill = "Fluorouracil", palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature relative contribution"))+
      scale_y_sqrt(labels = function(x) paste0(round(x*100), "%"),
                   limits = c(0,round(max(Y_pos))),
                   breaks=c(0.01,0.01,0.05,0.1,0.25,0.5,1))+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.02, colour="grey20", size=3.5)

  }
  
  
  if (mode=="absolute" & control_cohort=="No")
  {
    data <- data %>% filter(Fluorouracil_new != "Treated naive")
    data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Not 5-FU pretreated","5-FU pretreated"))
    data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
    data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)

    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo ~ Fluorouracil_new, data = data,method = "wilcox.test"));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- stats_data_plot[,c("V1")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    nn = data %>% group_by(Fluorouracil_new) %>% tally()
    
    set_max <- data %>% group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_log),median_sign = median(SIGN_5FU_denovo_log), SD_sign = sd(SIGN_5FU_denovo_log),max(SIGN_5FU_denovo_log), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_log)`)
    
    if(number<4){
      if(number==2){
        Y_pos=set_max+0.5
      }else{
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
      Y_pos = Y_pos+0.5
      }
    }
    
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "SIGN_5FU_denovo",
                      color = "black",fill = "Fluorouracil_new", palette = c( "#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      #stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      stat_compare_means(method = "wilcox.test", label.y = round(Y_pos),size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature absolute contribution"))+
      scale_y_continuous(trans='log10',
                         breaks=10**(1:round(max(Y_pos))),
                         limits = c(0.4,10^round(max(Y_pos))),
                         labels = scales::comma)+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.4, colour="grey20", size=4)
    
  }
  if (mode=="absolute" & control_cohort=="Yes")
  {
    data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
    data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)
    
    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo ~ Fluorouracil_new, data = data,method = "wilcox.test"));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    #stats_data_plot <- stats_data_plot[,c("V1")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    nn = data %>% group_by(Fluorouracil_new) %>% tally()
    
    set_max <- data %>% group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_log),median_sign = median(SIGN_5FU_denovo_log), SD_sign = sd(SIGN_5FU_denovo_log),max(SIGN_5FU_denovo_log), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_log)`)
    
    if(number<4){
      if(number==2){
        Y_pos=set_max+0.5
      }else{
        Y_pos = c(10/10,9/10,8/10)
        Y_pos = set_max*Y_pos
        Y_pos = Y_pos[1:number]
        Y_pos = Y_pos+0.75
      }
    }
    
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "SIGN_5FU_denovo",
                      color = "black",fill = "Fluorouracil_new", palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      #stat_compare_means(method = "wilcox.test", label.y = round(Y_pos),size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature absolute contribution"))+
      scale_y_continuous(trans='log10',
                         breaks=10**(1:round(max(Y_pos))),
                         limits = c(0.4,10^round(max(Y_pos))),
                         labels = scales::comma)+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.4, colour="grey20", size=4)
    
  }
  if(is.null(pdf.path)){
      return(plot)
    } else {
      pdf(pdf.path,5,5)
      plot(plot)
      dev.off()}
}
plot5Fucomparison(data = plotting_data_cohort,mode="relative") 
plot5Fucomparison(data = plotting_data_cohort,mode="absolute",control_cohort="No")
plot5Fucomparison(data = plotting_data_cohort,mode="absolute",control_cohort="Yes")
plot5Fucomparison(data = plotting_data_cohort,mode="relative",pdf.path=sprintf("%s%s_5-FU_relative_boxplot.pdf",dirpath,plot_name)) 
plot5Fucomparison(data = plotting_data_cohort,mode="absolute",control_cohort="No",pdf.path=sprintf("%s%s_5-FU_absolute_boxplot.pdf",dirpath,plot_name)) 
plot5Fucomparison(data = plotting_data_cohort,mode="absolute",control_cohort="Yes",pdf.path=sprintf("%s%s_5-FU_absolute_boxplot_with_control.pdf",dirpath,plot_name)) 

plot5Fucomparison_TGmuts <- function(data = plotting_data_cohort,mode="absolute", pdf.path = NULL){
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  #check number of samples
  
  
  
  if (mode=="relative")
  {
    #select_test
    #Kruskal???Wallis
    #res.dunn <- dunnTest(SIGN_5FU_denovo_rel ~ Fluorouracil, data = data, method = "holm")
    #control_compare <- which(grepl("0 -",res.dunn$res$Comparison))
    #res <- data.frame(res.dunn$res[control_compare,])
    
    #wilcox
    #wilcox <- wilcox.test(SIGN_5FU_denovo_rel ~ Fluorouracil, data = plotting_data_cohort[which(plotting_data_cohort$Cohort=="Metastatic"),])
    #print(wilcox$p.value)
    
    names(data)
    stats_data <- as.data.frame(compare_means(T2G_rel ~ Fluorouracil_new, data = data));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    nn = data %>% group_by(Fluorouracil_new) %>% tally()
    
    set_max <- data %>% group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(T2G_rel),median_sign = median(T2G_rel), SD_sign = sd(T2G_rel),max(T2G_rel), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- 2.5*max(set_max$`max(T2G_rel)`)
    set_max <- as.numeric(format(round(set_max, 1), nsmall = 1))
    
    if(number<4){
      Y_pos = c(98/100,89/100,84/100)
      Y_pos = set_max*Y_pos
      #Y_pos = 1*Y_pos
      Y_pos = Y_pos[1:number]
    }
    
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "T2G_rel",
                      color = "black",fill = "Fluorouracil", palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = 1.3*Y_pos,size = 4)+
      xlab(c(""))+
      ylab(c("relative contribution of \n C[T>G]T mutations"))+
      scale_y_sqrt(labels = function(x) paste0(round(x*100), "%"),
                   limits = c(0,set_max),
                   breaks=c(0.01,0.01,0.05,0.1,0.25,set_max))+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=0, colour="grey20", size=3.5)
    
  }
  
  
  if (mode=="absolute")
  {
    data <- data %>% filter(Fluorouracil_new != "Treated naive")
    data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Not 5-FU pretreated","5-FU pretreated"))
    data$T2G <- data$T2G+1
    data$T2G_log <- log10(data$SIGN_5FU_denovo)
    
    stats_data <- as.data.frame(compare_means(T2G ~ Fluorouracil_new, data = data,method = "wilcox.test"));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- stats_data_plot[,c("V1")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    nn = data %>% group_by(Fluorouracil_new) %>% tally()
    
    set_max <- data %>% group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(T2G_log),median_sign = median(T2G_log), SD_sign = sd(T2G_log),max(T2G_log), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(T2G_log)`)
    
    if(number<4){
      if(number==2){
        Y_pos=set_max
      }else{
        Y_pos = c(10/10,9/10,8/10)
        Y_pos = set_max*Y_pos
        Y_pos = Y_pos[1:number]
        Y_pos = Y_pos+0.5
      }
    }
    
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "T2G",
                      color = "black",fill = "Fluorouracil_new", palette = c( "#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      #stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      stat_compare_means(method = "wilcox.test", label.y = round(Y_pos),size = 4)+
      xlab(c(""))+
      ylab(c("absolute contribution of \n C[T>G]T mutations"))+
      scale_y_continuous(trans='log10',
                         breaks=10**(1:round(max(Y_pos))),
                         limits = c(0.4,10^round(max(Y_pos))),
                         labels = scales::comma)+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.4, colour="grey20", size=4)
    
  }
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()}
}
plot5Fucomparison_TGmuts(data = plotting_data_cohort,mode="relative",pdf.path=sprintf("%s%s_5-FU_TG_relative_boxplot.pdf",dirpath,plot_name)) 
plot5Fucomparison_TGmuts(data = plotting_data_cohort,mode="absolute",pdf.path=sprintf("%s%s_5-FU_TG_absolute_boxplot.pdf",dirpath,plot_name)) 


plottreatmenttypecomparison <- function(data = plotting_data_cohort,mode="relative", pdf.path = NULL){
  
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "No treatment", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","No treatment","5-FU pretreated"))
  names(data)
  #remove samples with multiple treatments
  data_subset=data[which(rowSums(as.matrix(data[c("Capecitabine","Fluoruracil","Folinic_acid")]),na.rm = T)<=1),]
  #remove control data
  data_subset <- data_subset %>% dplyr::filter(Cohort == "Metastatic") %>% dplyr::select(sampleId,Fluorouracil_new,Capecitabine,Fluoruracil,Folinic_acid)
  data_subset[is.na(data_subset)] <- 0
  data_m <- melt(data_subset)
  rm(data_subset)
  data_m <- data_m %>% dplyr::mutate(Treatment_new = ifelse(Fluorouracil_new == "No treatment", "No treatment", ifelse(Fluorouracil_new == "5-FU pretreated" &
                                                                                                                variable == "Capecitabine" &
                                                                                                                value == 1,"Capecitabine", ifelse(Fluorouracil_new == "5-FU pretreated" &
                                                                                                                                                    variable == "Fluoruracil" &
                                                                                                                                                    value == 1,"Fluoruracil", ifelse(Fluorouracil_new == "5-FU pretreated" &
                                                                                                                                                                                      variable == "Folinic_acid" &
                                                                                                                                                                                       value == 1,"Folinic_acid","to_remove")))))
  data_m$Treatment_new <- factor(data_m$Treatment_new, levels = c("No treatment","Fluoruracil","Capecitabine","Folinic_acid"))
  data_m <- data_m %>% dplyr::filter(Treatment_new != "to_remove")
  data_m <- data_m %>% dplyr::filter(Treatment_new != "Folinic_acid")
  data_m <- data_m %>% dplyr::select(sampleId,Treatment_new,value) %>% distinct()

  nn = data_m %>% dplyr::group_by(Treatment_new) %>% tally();nn
  
  data_m=dplyr::left_join(data_m,data[c("sampleId","SIGN_5FU_denovo_rel","SIGN_5FU_denovo")],by="sampleId")
  

  
  if (mode=="relative")
  {
    #select_test
    #Kruskal???Wallis
    #res.dunn <- dunnTest(SIGN_5FU_denovo_rel ~ Fluorouracil, data = data, method = "holm")
    #control_compare <- which(grepl("0 -",res.dunn$res$Comparison))
    #res <- data.frame(res.dunn$res[control_compare,])
    
    #wilcox
    #wilcox <- wilcox.test(SIGN_5FU_denovo_rel ~ Fluorouracil, data = plotting_data_cohort[which(plotting_data_cohort$Cohort=="Metastatic"),])
    #print(wilcox$p.value)
    
    
    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo_rel ~ Treatment_new, data = data_m));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    set_max <- data_m %>% dplyr::group_by(Treatment_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_rel),median_sign = median(SIGN_5FU_denovo_rel), SD_sign = sd(SIGN_5FU_denovo_rel),max(SIGN_5FU_denovo_rel), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_rel)`)
    
    if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
    }else{
      Y_pos = c(10/10,9.5/10,9/10,8.5/10,8/10,7.5/10,7/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
    }
    
    
    plot <- ggboxplot(data_m, x = "Treatment_new", y = "SIGN_5FU_denovo_rel",
                      color = "black",fill = "Treatment_new", palette = c("#E7B800","#FC4E07", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature relative contribution"))+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.015, colour="grey20", size=4)
    
  }
  
  
  if (mode=="absolute")
  {
    data_m$SIGN_5FU_denovo <- data_m$SIGN_5FU_denovo+1
    data_m$SIGN_5FU_denovo_log <- log10(data_m$SIGN_5FU_denovo)
    
    
    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo_log ~ Treatment_new, data = data_m));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    set_max <- data_m %>% dplyr::group_by(Treatment_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_log),median_sign = median(SIGN_5FU_denovo_log), SD_sign = sd(SIGN_5FU_denovo_log),max(SIGN_5FU_denovo_log), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_log)`)
    
    if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
      Y_pos = Y_pos+0.5
      
    }
    
    
    plot <- ggboxplot(data_m, x = "Treatment_new", y = "SIGN_5FU_denovo",
                      color = "black",fill = "Treatment_new", palette = c("#E7B800","#FC4E07", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature absolute contribution"))+
      scale_y_continuous(trans='log10',
                         breaks=10**(1:round(max(Y_pos))-1),
                         limits = c(0.4,10^round(max(Y_pos))),
                         labels = scales::comma)+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.4, colour="grey20", size=4)
    
  }
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()}
}
plottreatmenttypecomparison(data = plotting_data_cohort,mode="relative",pdf.path=sprintf("%s%s_treatmentType_relative_boxplot.pdf",dirpath,plot_name)) 
plottreatmenttypecomparison(data = plotting_data_cohort,mode="absolute",pdf.path=sprintf("%s%s_treatmentType_boxplot.pdf",dirpath,plot_name)) 

linear_fit_5FU_sign17 <- function(data = plotting_data_cohort, pdf.path = NULL){
  
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "No treatment", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","No treatment","5-FU pretreated"))
  
  fit <- lm(SIGN_17 ~ SIGN_5FU_denovo, data = data)
  plot <- ggplot(data, aes(x = SIGN_17, y = SIGN_5FU_denovo)) +
                 geom_point(aes(colour=Cohort),size=3) + 
                 stat_smooth(method = "lm", col = "#ca0020") +
                 scale_colour_manual(values = c("Metastatic"="#ef8a62",
                                                "Primary"="#67a9cf"),name="Cohort")+
                 labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                         "Intercept =",signif(fit$coef[[1]],5 ),
                         " Slope =",signif(fit$coef[[2]], 5),
                         " P =",signif(summary(fit)$coef[2,4], 5)))+
    xlab(c("absolute contribution COSMIC signature 17"))+
    ylab(c("absolute contribution 5-FU signature"))
    
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,7,5)
    plot(plot)
    dev.off()}
}
linear_fit_5FU_sign17(data = plotting_data_cohort,pdf.path=sprintf("%s%s_linear_fit_5FU_sign17.pdf",dirpath,plot_name)) 

linear_fit_5FU_T2G <- function(data = plotting_data_cohort, pdf.path = NULL){
  
  data <- data %>% dplyr::mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "No treatment", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","No treatment","5-FU pretreated"))
  
  data <- data %>% dplyr::filter(mut_load < 500000)
  
  fit <- lm(T2G ~ SIGN_5FU_denovo, data = data)
  plot <- ggplot(data, aes(x = T2G, y = SIGN_5FU_denovo)) +
    geom_point(aes(colour=Cohort),size=3) + 
    stat_smooth(method = "lm", col = "#ca0020") +
    scale_colour_manual(values = c("Metastatic"="#ef8a62",
                                   "Primary"="#67a9cf"),name="Cohort")+
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))+
    xlab(c("absolute contribution C[T>G]T"))+
    ylab(c("absolute contribution 5-FU signature"))
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,7,5)
    plot(plot)
    dev.off()}
}
linear_fit_5FU_T2G(data = plotting_data_cohort,pdf.path=sprintf("%s%s_linear_fit_5FU_sign17.pdf",dirpath,plot_name)) 


TYMS_treatment <- function(data = plotting_data_cohort, pdf.path = NULL){
  
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  names(data)
  #remove samples with multiple treatments
  #remove control data
  data <- data %>% dplyr::filter(Cohort == "Metastatic")
  nn = data %>% dplyr::group_by(Fluorouracil_new) %>% tally();nn
  
  print(data %>% dplyr::group_by(Fluorouracil_new) %>%
          dplyr::summarize(count = n(),mean_sign = mean(copyNumber_TYMS),median_sign = median(copyNumber_TYMS), SD_sign = sd(copyNumber_TYMS),max(copyNumber_TYMS), sum(median_sign,SD_sign)) ) 
  
  print(data %>% dplyr::filter(copyNumber_TYMS>4) %>% dplyr::group_by(Fluorouracil_new) %>%
          dplyr::summarize(count = n(),mean_sign = mean(copyNumber_TYMS),median_sign = median(copyNumber_TYMS), SD_sign = sd(copyNumber_TYMS),max(copyNumber_TYMS), sum(median_sign,SD_sign)) )
  
   
  #fisher on drivers
  #treated or not
  pathos_nonFU <- as.numeric(8)
  samples_nonFU <- as.numeric(121)
  pathos_FU <- as.numeric(44)
  samples_FU <-  as.numeric(231)
  
  challenge.df = matrix(c(pathos_FU,pathos_nonFU,samples_FU,samples_nonFU), nrow = 2)
  fisher.test(challenge.df, alternative ="greater")
  
  stats_data <- as.data.frame(compare_means(copyNumber_TYMS ~ Fluorouracil_new, data = data));print(stats_data)
  stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
  #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
  stats_data_plot <- lapply(stats_data_plot, as.character)
  stats_data_plot <- as.list(stats_data_plot)
  
  set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
    dplyr::summarize(count = n(),mean_sign = mean(copyNumber_TYMS),median_sign = median(copyNumber_TYMS), SD_sign = sd(copyNumber_TYMS),max(copyNumber_TYMS), sum(median_sign,SD_sign)) 
  set_max <- as.data.frame(set_max)
  number  <- nrow(set_max)
  set_max <- max(set_max$`max(copyNumber_TYMS)`)
  
  if(number<4){
    Y_pos = c(10/10,9/10,8/10)
    Y_pos = set_max*Y_pos
    Y_pos = Y_pos[1:number]
    Y_pos = Y_pos+1
    
  }
  
  plot <- ggboxplot(data, x = "Fluorouracil_new", y = "copyNumber_TYMS",
                    color = "black",fill = "Fluorouracil_new", palette = c("#008837","#7b3294"))+ 
    #stat_compare_means(aes(label = ..p.signif..))+
    stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos+1,size = 4)+
    xlab(c(""))+
    ylab(c("TYMS copy number level"))+
    scale_y_continuous(limits = c(-2,round(max(Y_pos))+5),labels = scales::comma)+
    #scale_y_continuous(limits = c(0, 0.5))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")+
    geom_text(data=nn, aes(label=paste0("n=", nn$n)),y=-1 ,colour="grey20", size=4)

  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,4,5)
    plot(plot)
    dev.off()}
}
TYMS_treatment(data = plotting_data_cohort,pdf.path=sprintf("%s%s_TYMS_treatment.pdf",dirpath,plot_name)) 

TYMS_signature <- function(data = plotting_data_cohort, pdf.path = NULL){
  
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  names(data)
  #remove samples with multiple treatments
  #remove control data
  data <- data %>% filter(Cohort == "Metastatic")
  #remove not_pretreated data
  data <- data %>% filter(Fluorouracil_new == "5-FU pretreated")
  data$SIGN_5FU_denovo_normalized <- data$SIGN_5FU_denovo / data$mut_load
  
  #exponential.model <- lm(log10(SIGN_5FU_denovo+1)~ copyNumber_TYMS, data = data)
  #summary(exponential.model)
  #copyNumber_TYMS_model <- seq(0, 25, 0.1)
  #
  #newframe <- data.frame(copyNumber_TYMS = seq(min(data$copyNumber_TYMS), max(data$copyNumber_TYMS), length = 1000))
  #pred <- predict(exponential.model, newdata = newframe, interval = "prediction")
  #pred <- cbind(pred,newframe)
#
  #ggplot()+
  #  geom_point(data = data,
  #             aes(x = copyNumber_TYMS,
  #                 y =  SIGN_5FU_denovo),
  #             size=3,color="#7b3294")+
  #  geom_line(data = pred,
  #            aes(x = copyNumber_TYMS,
  #                y =  exp(fit)),color="#7b3294")

  plot <- ggplot(data = data,
         aes(x = copyNumber_TYMS,
             y =  SIGN_5FU_denovo))+
    geom_point(size=3,color="#7b3294",alpha = 0.7)+ 
    xlab(c("TYMS copy number level"))+
    ylab(c("5-FU signature absolute contribution"))+
    scale_x_continuous(breaks = c(0,1,2,3,4,5,10,15,20,25,30),labels = scales::comma)+
    #scale_y_continuous(limits = c(0, 0.5))+
    theme(legend.position = "none")

  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,7,5)
    plot(plot)
    dev.off()}
}
TYMS_signature(data = plotting_data_cohort,pdf.path=sprintf("%s%s_TYMS_signature.pdf",dirpath,plot_name)) 

plot5Fu_TMB_comparison <- function(data = plotting_data_cohort,control_cohort="No",TMB="Yes",pdf.path = NULL){
  data <- data %>% dplyr::mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  if (control_cohort=="No"&TMB=="No"){
  #remove control data
  data <- data %>% dplyr::filter(Cohort == "Metastatic")
  #check number of samples
  nn = data %>% dplyr::group_by(Fluorouracil_new) %>% tally()
  
  data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
  data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)
  data$mut_load <- data$mut_load+1
  data$mut_load_log <- log10(data$mut_load)
    
  stats_data <- as.data.frame(compare_means(mut_load_log ~ Fluorouracil_new, data = data));print(stats_data)
  stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
  #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
  stats_data_plot <- lapply(stats_data_plot, as.character)
  stats_data_plot <- as.list(stats_data_plot)
  
  print(set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
          dplyr::summarize(count = n(),mean_sign = mean(mut_load),median_sign = median(mut_load), SD_sign = sd(mut_load),max(mut_load), sum(median_sign,SD_sign)) )
  
  set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
    dplyr::summarize(count = n(),mean_sign = mean(mut_load_log),median_sign = median(mut_load_log), SD_sign = sd(mut_load_log),max(mut_load_log), sum(median_sign,SD_sign)) 
  set_max <- as.data.frame(set_max)
  number  <- nrow(set_max)
  set_max <- max(set_max$`max(mut_load_log)`)
    
  if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
      Y_pos = Y_pos+1
      
    }
    max_new <- round(max(data$mut_load)/10000)*10000
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "mut_load",
                      color = "black",fill = "Fluorouracil_new", palette = c("#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = max_new,size = 4)+
      xlab(c(""))+
      ylab(c("Tumor mutation burden"))+
      scale_y_continuous(limits = c(0,max_new),labels = scales::comma)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=3, colour="grey20", size=4)
  }
  if (control_cohort=="No"&TMB=="Yes"){
    #remove control data
    data <- data %>% dplyr::filter(Cohort == "Metastatic")
    #check number of samples
    nn = data %>% dplyr::group_by(Fluorouracil_new) %>% tally()

    
    stats_data <- as.data.frame(compare_means(TMB ~ Fluorouracil_new, data = data));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    print(set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
            dplyr::summarize(count = n(),mean_sign = mean(TMB),median_sign = median(TMB), SD_sign = sd(TMB),max(TMB), sum(median_sign,SD_sign)) )
    
    set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(TMB),median_sign = median(TMB), SD_sign = sd(TMB),max(TMB), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(TMB)`)
    
    if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
      Y_pos = Y_pos+1
      
    }
    max_new <- 11
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "TMB",
                      color = "black",fill = "Fluorouracil_new", palette = c("#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = 11,size = 4)+
      xlab(c(""))+
      ylab(c("Tumor mutation burden"))+
      scale_y_continuous(limits = c(-0.8,12),
                         breaks=c(0,2.5,5,7.5,10))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.5, colour="grey20", size=4)
  }
  if (control_cohort=="Yes"){
    #check number of samples
    nn = data %>% dplyr::group_by(Fluorouracil_new) %>% tally()
    
    data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
    data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)
    data$mut_load <- data$mut_load+1
    data$mut_load_log <- log10(data$mut_load)
    
    stats_data <- as.data.frame(compare_means(mut_load_log ~ Fluorouracil_new, data = data));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    print(set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
            dplyr::summarize(count = n(),mean_sign = mean(mut_load),median_sign = median(mut_load), SD_sign = sd(mut_load),max(mut_load), sum(median_sign,SD_sign)) )
    
    set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
      dplyr::summarize(count = n(),mean_sign = mean(mut_load_log),median_sign = median(mut_load_log), SD_sign = sd(mut_load_log),max(mut_load_log), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(mut_load_log)`)
    
    if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
      Y_pos = Y_pos+1
      
    }
    #max_new <- round(max(data$mut_load)/10000)*10000
    plot <- ggboxplot(data, x = "Fluorouracil_new", y = "mut_load",
              color = "black",fill = "Fluorouracil_new", palette = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      #stat_compare_means(method = "wilcox.test", label.y = round(Y_pos),size = 4)+
      xlab(c(""))+
      ylab(c("Tumor mutation burden"))+
      scale_y_continuous(trans='log10',
                         breaks=10**(1:round(max(Y_pos))),
                         limits = c(100,10^round(max(Y_pos+0.1))),
                         labels = scales::comma)+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=2, colour="grey20", size=4)
  }

  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()}
}
plot5Fu_TMB_comparison(data = plotting_data_cohort,control_cohort="No") 
plot5Fu_TMB_comparison(data = plotting_data_cohort,control_cohort="Yes") 
plot5Fu_TMB_comparison(data = plotting_data_cohort,control_cohort="No",pdf.path=sprintf("%s%s_plot5Fu_TMB_comparison_WO_control.pdf",dirpath,plot_name)) 
plot5Fu_TMB_comparison(data = plotting_data_cohort,control_cohort="Yes",pdf.path=sprintf("%s%s_plot5Fu_TMB_comparison_W_control.pdf",dirpath,plot_name)) 


plot5Fu_mut_load_contribution_per_patient <- function(data = plotting_data_cohort,binned="No",pdf.path = NULL){
  data <- data %>% dplyr::mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  #remove control data
  bins<-c("0% - 10%","10% - 20%","20% - 30%","30% - 40%","40% - 50%","50% - 100%")
  data <- data %>% dplyr::filter(Cohort == "Metastatic", Fluorouracil_new=="5-FU pretreated") %>%
    dplyr::mutate(bin = ifelse(SIGN_5FU_denovo_rel <= 0.1, "0% - 10%", 
                               ifelse(SIGN_5FU_denovo_rel <= 0.2, "10% - 20%",
                                      ifelse(SIGN_5FU_denovo_rel <= 0.3, "20% - 30%",
                                             ifelse(SIGN_5FU_denovo_rel <= 0.4, "30% - 40%",
                                                    ifelse(SIGN_5FU_denovo_rel <= 0.5, "40% - 50%","50% - 100%"))))))
  hist(data$mut_load)
  data <- data %>%
    dplyr::mutate(mut_load_bin = ifelse(mut_load <= 10000, "0 - 10", 
                               ifelse(mut_load <= 20000, "10 - 20",
                                      ifelse(mut_load <= 30000, "20 - 30","> 30"))))
  max(data$SIGN_5FU_denovo_rel)
  data$bin <- factor(data$bin, levels = bins)
  #check number of samples
  nn = data %>% dplyr::group_by(bin) %>% tally() %>% as.data.frame()
  nn2 = setNames(data.frame(matrix(ncol = 1, nrow = length(bins))),c("bin"))
  nn2$bin <- bins
  nn <- dplyr::full_join(nn,nn2, by="bin") %>% dplyr::mutate(n = ifelse(is.na(n), 0, n))
  print(nn)
  
  plot_height=max(nn$n)+10
  text_height=plot_height-2
  if (binned=="Yes"){
  plot <- ggplot(data, aes(x = bin)) +
    #geom_bar(fill = c("#FC4E07")) + 
    geom_bar(aes(fill = mut_load_bin)) + 
    ylab("5-FU contribution on TMB") + 
    scale_y_continuous(limits = c(0,plot_height),labels = scales::comma)+
    scale_x_discrete(drop = FALSE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),
          axis.text.y.left = element_text(size = 10),
          axis.ticks.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_rect(fill=NA),
          legend.position = "bottom")+
    geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=text_height, colour="black", size=4)
  }
  if (binned=="No"){
  plot <- ggplot(data, aes(x = bin)) +
    geom_bar(fill = c("#FC4E07")) + 
    #geom_bar(aes(fill = mut_load_bin)) + 
    ylab("5-FU contribution on TMB") + 
    scale_y_continuous(limits = c(0,plot_height),labels = scales::comma)+
    scale_x_discrete(drop = FALSE)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),
          axis.text.y.left = element_text(size = 10),
          axis.ticks.x=element_blank(),
          panel.background=element_blank(),
          panel.border=element_rect(fill=NA),
          legend.position = "bottom")+
    geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=text_height, colour="black", size=4)
  }
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()
    }
}
plot5Fu_mut_load_contribution_per_patient(data = plotting_data_cohort,binned="No") 
plot5Fu_mut_load_contribution_per_patient(data = plotting_data_cohort,binned="Yes",pdf.path=sprintf("%s%s_plot5Fu_mut_load_contribution_per_patient_binned.pdf",dirpath,plot_name)) 
plot5Fu_mut_load_contribution_per_patient(data = plotting_data_cohort,binned="No",pdf.path=sprintf("%s%s_plot5Fu_mut_load_contribution_per_patient.pdf",dirpath,plot_name)) 


plot_TMB_wo5FU_comparison <- function(data = plotting_data_cohort,pdf.path = NULL){
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  #remove control data
  data <- data %>% dplyr::filter(Cohort == "Metastatic")
  
  #check number of samples
  nn = data %>% dplyr::group_by(Fluorouracil_new) %>% tally()
  
  data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
  data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)
  data$mut_load <- data$mut_load+1
  data$mut_load_log <- log10(data$mut_load)
  
  data$mut_load_new <- data$mut_load - data$SIGN_5FU_denovo
  data$mut_load_new <- data$mut_load_new+1
  data$mut_load_new_log <- log10(data$mut_load_new)
  
  stats_data <- as.data.frame(compare_means(mut_load_new_log ~ Fluorouracil_new, data = data));print(stats_data)
  stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
  #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
  stats_data_plot <- lapply(stats_data_plot, as.character)
  stats_data_plot <- as.list(stats_data_plot)
  
  print(data %>% dplyr::group_by(Fluorouracil_new) %>%
          dplyr::summarize(count = n(),mean_sign = mean(mut_load_new),median_sign = median(mut_load_new), SD_sign = sd(mut_load_new),max(mut_load_new), sum(median_sign,SD_sign)) )
  set_max <- data %>% dplyr::group_by(Fluorouracil_new) %>%
    dplyr::summarize(count = n(),mean_sign = mean(mut_load_new_log),median_sign = median(mut_load_new_log), SD_sign = sd(mut_load_new_log),max(mut_load_new_log), sum(median_sign,SD_sign)) 
  set_max <- as.data.frame(set_max)
  number  <- nrow(set_max)
  set_max <- max(set_max$`max(mut_load_new_log)`)
  
  if(number<4){
    Y_pos = c(10/10,9/10,8/10)
    Y_pos = set_max*Y_pos
    Y_pos = Y_pos[1:number]
    Y_pos = Y_pos+1
    
  }
  
  max_new <- round(max(data$mut_load_new)/10000)*10000
  plot <- ggboxplot(data, x = "Fluorouracil_new", y = "mut_load_new",
                    color = "black",fill = "Fluorouracil_new", palette = c("#E7B800", "#FC4E07"))+ 
    #stat_compare_means(aes(label = ..p.signif..))+
    stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = max_new,size = 4)+
    xlab(c(""))+
    ylab(c("Tumor mutation burden"))+
    scale_y_continuous(limits = c(0,max_new),labels = scales::comma)+
    #scale_y_continuous(limits = c(0, 0.5))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")+
    geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=3, colour="grey20", size=4)
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()}
}
#plot_TMB_wo5FU_comparison(data = plotting_data_cohort,pdf.path=sprintf("%s%s_plot_TMB_WO_5FU_comparison.pdf",dirpath,plot_name)) 

plot_TMB_With_and_WO_5FU_comparison <- function(data = plotting_data_cohort,pdf.path = NULL){
  data <- data %>% dplyr::mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))
  #remove control data
  data <- data %>% dplyr::filter(Cohort == "Metastatic")
  #data <- data %>% dplyr::filter(hasSystemicPreTreatment=="Yes")

  data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
  data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)
  
  data$TMB_with_NMF_5FU <- data$mut_load+1
  data$mut_load_log <- log10(data$TMB_with_NMF_5FU)
  
  data$TMB_without_NMF_5FU <- data$TMB_with_NMF_5FU - data$SIGN_5FU_denovo
  data$TMB_without_NMF_5FU <- data$TMB_without_NMF_5FU+2
  data$mut_load_new_log <- log10(data$TMB_without_NMF_5FU)
  
  data <- data %>% dplyr::select(Fluorouracil_new,TMB_without_NMF_5FU,TMB_with_NMF_5FU,sampleId)
  data <- melt(data)
  data$value <- round(as.numeric(data$value))
  
  print(data %>% dplyr::group_by(Fluorouracil_new,variable) %>%
          dplyr::summarize(count = n(),mean_sign = mean(value),median_sign = median(value), SD_sign = sd(value),max(value), sum(median_sign,SD_sign)) )
  set_max <- data %>% dplyr::group_by(Fluorouracil_new,variable) %>%
    dplyr::summarize(count = n(),mean_sign = mean(value),median_sign = median(value), SD_sign = sd(value),max(value), sum(median_sign,SD_sign)) 
  set_max <- as.data.frame(set_max)
  stats_data_plot <- as.data.frame(t(set_max[,c("Fluorouracil_new", "variable")]))
  stats_data_plot <- lapply(stats_data_plot, as.character)
  stats_data_plot <- as.list(stats_data_plot)
  number  <- nrow(set_max)
  set_max <- max(set_max$`max(value)`)
  set_max <- log10(set_max)
  
  #check number of samples
  nn = data %>% dplyr::group_by(Fluorouracil_new,variable) %>% tally()
  
  if(number<4){
    Y_pos = c(10/10,9/10,8/10)
    Y_pos = set_max*Y_pos
    Y_pos = Y_pos[1:number]
    Y_pos = Y_pos+1
  }else{
    Y_pos = c(10/10,9.5/10,9/10,8.5/10,8/10,7.5/10,7/10)
    Y_pos = set_max*Y_pos
    Y_pos = Y_pos[1:number]
  }
  
  data$value_log <- log10(data$value)
  
  fit1.lme <- lme(value~variable*Fluorouracil_new,data = data,random=~+1|sampleId)
  anova(fit1.lme )
  fitlist <- as.data.frame(anova(lme(value~variable*Fluorouracil_new,data = data,random=~+1|sampleId)))$`p-value`
  pvalue_linear_mixed_model <- tail(fitlist, n=1)
  pvalue_linear_mixed_model <- format(pvalue_linear_mixed_model, scientific = T, digits = 2)
  pvalue_linear_mixed_model <- paste0("ANOVA linear mixed model; P=",pvalue_linear_mixed_model)
  
  

  plot <- ggboxplot(data, x = "variable", y = "value",
                    color = "black",facet.by = "Fluorouracil_new",fill="variable",
                    title=pvalue_linear_mixed_model) + 
    #stat_compare_means(aes(label = ..p.signif..))+
    stat_compare_means(method = "wilcox.test", size = 4)+
    xlab(c(""))+
    #ylab( expression(log[10]~(mutational~load)))+
    ylab("Mutation contribution")+
    scale_y_continuous(limits = c(0,round(max(data$value)/10000)*10000),labels = scales::comma)+
    #scale_y_continuous(limits = c(3, 5))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")+
    geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=3, colour="grey20", size=4)
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()}
}
plot_TMB_With_and_WO_5FU_comparison(data = plotting_data_cohort) 
plot_TMB_With_and_WO_5FU_comparison(data = plotting_data_cohort,pdf.path=sprintf("%s%s_plot_TMB_With_and_WO_5FU_comparison.pdf",dirpath,plot_name)) 

################
################
##SUBCLONAL mutations
################
################
load("~/surfdrive/Shared/Sig17/HMF_data/cohort_analyse/clonality_table_colon_breast_NMF16.RData")




clonality_table <- clonality_table
clonality_table <- clonality_table_colon
clonality_table <- clonality_table_breast
names(clonality_table)
clonality_table <- clonality_table[which(clonality_table$sampleId %in% plotting_data_cohort$sampleId),]
sampleIDs_with_subclonal_muts <- clonality_table %>% dplyr::filter(clonality=="SUBCLONAL"& mut_load>500) %>% dplyr::pull(sampleId)
clonality_table <- clonality_table[which(clonality_table$sampleId %in% sampleIDs_with_subclonal_muts),]
clonality_table=dplyr::inner_join(clonality_table,plotting_data_cohort[which(plotting_data_cohort$Cohort=="Metastatic"),][c("sampleId","Fluorouracil")],by="sampleId")
length(sampleIDs_with_subclonal_muts)
clonality_table_transformed=data.frame(sampleId=clonality_table$sampleId,
                                       Clonality=clonality_table$clonality,
                                       SIGN_5FU_denovo=clonality_table$NMF_H,    #breast: NMF_breast_F    colon: NMF_colon_G
                                       SIGN_5FU_denovo_rel=clonality_table$NMF_H_rel,
                                       SIGN_17=clonality_table$Signature.17,
                                       SIGN_17_rel=clonality_table$Signature.17_rel,
                                       T2G=clonality_table$`C[T>G]T`,
                                       T2G_rel=clonality_table$`C[T>G]T_rel`,
                                       Fluorouracil=as.factor(clonality_table$Fluorouracil),
                                       mut_load=clonality_table$mut_load)

plot_clonal_muts <- function(data = clonality_table_transformed,mode="relative", pdf.path = NULL){
  data <- data %>% dplyr::mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Not 5-FU pretreated","5-FU pretreated"))
  #check number of samples
  data <- data %>% dplyr::filter(Fluorouracil_new=="5-FU pretreated")
  nn = data  %>% dplyr::group_by(Fluorouracil_new,Clonality) %>% tally();nn
  
  
  if (mode=="relative")
  {
    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo_rel ~ Clonality, data = data));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    set_max <- data %>% dplyr::group_by(Clonality) %>%
      dplyr::summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_rel),median_sign = median(SIGN_5FU_denovo_rel), SD_sign = sd(SIGN_5FU_denovo_rel),max(SIGN_5FU_denovo_rel), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_rel)`)
    
    if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
    }
    
    
    plot <- ggboxplot(data, x = "Clonality", y = "SIGN_5FU_denovo_rel",
                      color = "black",fill = "Clonality", palette = c("#a6611a", "#7570b3"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = 0.98,size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature relative contribution"))+
      scale_y_sqrt(labels = function(x) paste0(round(x*100), "%"),
                   limits = c(0,round(max(Y_pos))),
                   breaks=c(0.01,0.01,0.05,0.1,0.25,0.5,1))+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.019, colour="grey20", size=4)
    
  }
  
  
  if (mode=="absolute")
  {
    data$SIGN_5FU_denovo <- data$SIGN_5FU_denovo+1
    data$SIGN_5FU_denovo_log <- log10(data$SIGN_5FU_denovo)
    
    stats_data <- as.data.frame(compare_means(SIGN_5FU_denovo_log ~ Clonality, data = data));print(stats_data)
    stats_data_plot <- as.data.frame(t(stats_data[,c("group1", "group2")]))
    #stats_data_plot <- stats_data_plot[,c("V2","V1","V3")]
    stats_data_plot <- lapply(stats_data_plot, as.character)
    stats_data_plot <- as.list(stats_data_plot)
    
    set_max <- data %>% group_by(Clonality) %>%
      summarize(count = n(),mean_sign = mean(SIGN_5FU_denovo_log),median_sign = median(SIGN_5FU_denovo_log), SD_sign = sd(SIGN_5FU_denovo_log),max(SIGN_5FU_denovo_log), sum(median_sign,SD_sign)) 
    set_max <- as.data.frame(set_max)
    number  <- nrow(set_max)
    set_max <- max(set_max$`max(SIGN_5FU_denovo_log)`)
    
    if(number<4){
      Y_pos = c(10/10,9/10,8/10)
      Y_pos = set_max*Y_pos
      Y_pos = Y_pos[1:number]
      Y_pos = Y_pos+1
      
    }
    
    
    plot <- ggboxplot(data, x = "Clonality", y = "SIGN_5FU_denovo",
                      color = "black",fill = "Clonality", palette = c("#a6611a", "#7570b3"))+ 
      #stat_compare_means(aes(label = ..p.signif..))+
      stat_compare_means(method = "wilcox.test",comparisons = stats_data_plot, label.y = Y_pos,size = 4)+
      xlab(c(""))+
      ylab(c("5-FU signature absolute contribution"))+
      scale_y_continuous(trans='log10',
                         breaks=10**(1:round(max(Y_pos))-1),
                         limits = c(0.4,10^round(max(Y_pos))),
                         labels = scales::comma)+
      #scale_y_continuous(limits = c(0, 0.5))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")+
      geom_text(data=nn, aes(label=paste0("n=", nn$n)),  y=-0.4, colour="grey20", size=4)
    
  }
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,5,5)
    plot(plot)
    dev.off()}
}
plot_clonal_muts(data = clonality_table_transformed,mode="relative",pdf.path=sprintf("%s%s_clonality_relative_boxplot.pdf",dirpath,plot_name)) 

CheckTMB <- function(data = somatic_clinical_subset, pdf.path = NULL){
  data_meta <- data %>% filter(Cohort == "Metastatic",hasSystemicPreTreatment=="Yes",TMB <=10) %>% dplyr::select(Fluorouracil,NMF_A_rel:NMF_P_rel) #NMF_colon_A:NMF_colon_J  #NMF_breast_A:NMF_breast_J
  data_prim <- data %>% filter(Cohort == "Primary",TMB <=10) %>% dplyr::select(Fluorouracil,NMF_A_rel:NMF_P_rel) #NMF_colon_A:NMF_colon_J  #NMF_breast_A:NMF_breast_J

  data <- rbind(data_meta,data_prim)
  data <- data %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  
  data$Fluorouracil_new <- factor(data$Fluorouracil_new, levels = c("Treated naive","Not 5-FU pretreated","5-FU pretreated"))


  data <- data %>% dplyr::select(Fluorouracil_new,NMF_A_rel:NMF_P_rel)
  data2 = data

  wilcoxresults <- as.list(lapply(names(data)[-1],function(x)
    as.data.frame(compare_means(as.formula(paste(x,"Fluorouracil_new",sep="~")),data=data2,method = "wilcox.test"))))
  wilcox.test.pval <- sapply(wilcoxresults, '[[', 'p')
  print(rbind(names(data)[-1],wilcox.test.pval))
  
  checkTMB_m = melt(data,id=c("Fluorouracil"))
  checkTMB_m$sign_contribution=log10(checkTMB_m$value+1)
  checkTMB_m <- checkTMB_m %>% mutate(Fluorouracil_new = ifelse(Fluorouracil == "0", "Not 5-FU pretreated", ifelse(Fluorouracil == "1", "5-FU pretreated", "Treated naive")))
  checkTMB_m$Fluorouracil_new <- factor(checkTMB_m$Fluorouracil_new, levels = c("Not 5-FU pretreated","5-FU pretreated"))
  
  plot <- ggboxplot(checkTMB_m, x = "Fluorouracil_new", y = "sign_contribution",
                    color = "black",fill = "Fluorouracil_new",facet.by="variable")+  #fill = "variable", #palette = c("#00AFBB", "#E7B800", "#FC4E07")
    #stat_compare_means(aes(label = ..p.signif..))+
    stat_compare_means(method = "wilcox.test")+
    xlab(c(""))+
    ylab(c("Signature mutation contribution"))+
    #scale_y_continuous(limits = c(0, 0.5))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none")
  
  if(is.null(pdf.path)){
    return(plot)
  } else {
    pdf(pdf.path,8,10)
    plot(plot)
    dev.off()
  }
  
}

