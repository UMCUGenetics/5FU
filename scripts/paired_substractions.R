library(SPARQL, quietly = T)
library(stringr, quietly = T)
library(Biostrings, quietly = T)
library(dplyr,quietly = T)
library(foreach, quietly = T)

options(stringsAsFactors = F)

endpoint <- "http://localhost:8890/sparql-auth"
auth_options <- curlOptions(userpwd="avanhoeck:XXXXXXXXX")

connect <- function(){
  system2(command="ssh", 
          args=c("-f",
                 "-o",
                 "ExitOnForwardFailure=yes",
                 "-i",
                 "/home/cog/bvanderroest/.ssh/fedor13key",
                 "-L",
                 "8890:localhost:8890",
                 "fedor13",
                 "sleep",
                 "10"))
}

#========= Query =========#
args = commandArgs(trailingOnly=TRUE)

sprintf("START PROGRAM FOR %s", sampleid)

query <- sprintf(
  "
  PREFIX col: <http://sparqling-genomics/table2rdf/Column/>
  
  PREFIX sampleid_pf: <http://sparqling-genomics/Sample/>
  PREFIX chromosome_pf: <http://rdf.biosemantics.org/data/genomeassemblies/hg19>
  PREFIX clonality_pf: <http://sparqling-genomics/Clonality/>
  PREFIX filter_pf: <http://sparqling-genomics/FilterItem/>
  PREFIX type_pf: <http://sparqling-genomics/VariantCallType/>
  
  SELECT 
  STRAFTER(STR(?sampleid), \"http://sparqling-genomics/Sample/\") AS ?sampleid
  STRAFTER(STR(?chromosome), \"http://rdf.biosemantics.org/data/genomeassemblies/hg19#\") AS ?chromosome
  ?position
  ?ref
  ?alt
  STRAFTER(STR(?clonality), \"http://sparqling-genomics/Clonality/\") AS ?clonality
  ?trinucleotidecontext
  STRAFTER(STR(?type), \"http://sparqling-genomics/VariantCallType/\") AS ?type
  STRAFTER(STR(?filter), \"http://sparqling-genomics/FilterItem/\") AS ?filter
  
  FROM <http://hmfpatients/somaticvariant>
  WHERE {
  ?row col:sampleid ?sampleid .
  ?row col:chromosome ?chromosome .
  ?row col:position ?position .
  ?row col:ref ?ref .
  ?row col:alt ?alt .
  ?row col:clonality ?clonality .
  ?row col:trinucleotidecontext ?trinucleotidecontext .
  ?row col:filter ?filter .
  ?row col:type ?type .
  
  FILTER (?filter = filter_pf:PASS)
  FILTER (?type IN (type_pf:SNP))
  FILTER (?sampleid IN (sampleid_pf:%s))
  }
  
  ",
  sampleid)

connect()
query_out <- SPARQL(url = endpoint, curl_args = auth_options, query = query)$results

sprintf("QUERY DONE FOR %s", sampleid)

sample_name <- sampleid

if (endsWith(sample_name, "T")){
  system2(command = "mkdir", args = c(sprintf("/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/analysis/paired_samples/%s",sample_name)))
  sample_dir <- sample_name
} else if (endsWith(sample_name, "I")){
  split <- unlist(strsplit(sample_name, "TI"))[1]
  sample_dir <- paste0(split, "T")
  system2(command = "mkdir", args = c(sprintf("/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/analysis/paired_samples/%s",sample_dir)))
  
}

write.table(query_out, 
            sprintf("/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/analysis/paired_samples/%s/%s_ALL.txt",
                    sample_dir, sample_name),
            row.names = F)
write.table(subset(query_out, clonality == "CLONAL"),
            sprintf("/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/analysis/paired_samples/%s/%s_CLONAL.txt",
                    sample_dir, sample_name),
            row.names = F)
write.table(subset(query_out, clonality == "SUBCLONAL"),
            sprintf("/hpc/cog_bioinf/cuppen/project_data/HMF_data/DR010-update/analysis/paired_samples/%s/%s_SUBCLONAL.txt",
                    sample_dir, sample_name),
            row.names = F)

sprintf("PROGAM DONE FOR %s", sampleid)