loadMutationMatrix <- function(path, pattern, clonality = "all"){
    path_c = paste0(path,"snv_HMF_update_clonal/SNV_CLONAL")
    path_sc = paste0(path,"snv_HMF_update_subclonal/SNV_SUBCLONAL")
    
    filenames_c=list.files(path=path_c, pattern = pattern, full.names=TRUE)
    filenames_sc=list.files(path=path_sc, pattern = pattern, full.names=TRUE)
    
    datalist = lapply(filenames_c, function(x){read.table(file=x, header = T, sep="", row.names = 1)})
    datalist_HMF_data = do.call(cbind, datalist)
    mut_matrix_all=cbind(datalist_HMF_data)
    select_all <- c(df_final$sampleID)
    select_int = intersect(select_all,colnames(mut_matrix_all))
    mut_matrix_all <- mut_matrix_all[,select_int]
    ncol(mut_matrix_all)
    mut_matrix_clonal <- mut_matrix_all
    
    datalist = lapply(filenames_sc, function(x){read.table(file=x, header = T, sep="", row.names = 1)})
    datalist_HMF_data = do.call(cbind, datalist)
    mut_matrix_all=cbind(datalist_HMF_data)
    select_all <- c(df_final$sampleID)
    select_int = intersect(select_all,colnames(mut_matrix_all))
    mut_matrix_all <- mut_matrix_all[,select_int]
    ncol(mut_matrix_all)
    mut_matrix_subclonal <- mut_matrix_all
    
    high_sub_clonal <- colnames(mut_matrix_subclonal)[which(colSums(mut_matrix_subclonal) > 200 )]
    
    if (clonality == "CLONAL"){
      return(mut_matrix_clonal[,high_sub_clonal])
    } else if (clonality == "SUBCLONAL"){
      return(mut_matrix_subclonal[,high_sub_clonal])
    } else if (clonality == "all"){
      return(rbind(mut_matrix_clonal[,high_sub_clonal], mut_matrix_subclonal[,high_sub_clonal]))
    }
  
}
