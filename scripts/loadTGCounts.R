loadTGCounts <- function(path, pattern, df_treatment, clonality){
  
  mut_matrix <- loadMutationMatrix(path=path,
                                   pattern=pattern,
                                   clonality = "all")
  
  TG_counts_clonal <- mut_matrix["C[T>G]T",]
  TG_counts_clonal <- TG_counts_clonal / colSums(mut_matrix)
  TG_counts_clonal <- melt(TG_counts_clonal)
  TG_counts_clonal$clonality <- "CLONAL"
  
  TG_counts_subclonal <- mut_matrix["C[T>G]T1",]
  TG_counts_subclonal <- TG_counts_subclonal / colSums(mut_matrix)
  TG_counts_subclonal <- melt(TG_counts_subclonal)
  TG_counts_subclonal$clonality <- "SUBCLONAL"
  
  if (clonality == "all"){
    TG_counts <- rbind(TG_counts_clonal, TG_counts_subclonal)
  } else if (clonality == "CLONAL"){
    TG_counts <- TG_counts_clonal
  } else if (clonality == "SUBCLONAL"){
    TG_counts <- TG_counts_subclonal
  }
  
  TG_counts$mutation <- "T>G"
  TG_counts$variable <- as.vector(TG_counts$variable)
  TG_counts <- TG_counts[,c(4,1,2,3)]
  
}
