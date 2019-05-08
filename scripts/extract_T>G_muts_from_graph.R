library(SPARQL, quietly = T)
library(stringr, quietly = T)
library(Biostrings, quietly = T)
library(dplyr,quietly = T)
library(foreach, quietly = T)

options(stringsAsFactors = F)

endpoint <- "http://localhost:8890/sparql-auth"
auth_options <- curlOptions(userpwd="avanhoeck:XXXXXXXXXX")


#========= Query =========#



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
  STRAFTER(STR(?ref), \"http://sparqling-genomics/Sequence/\") AS ?ref
  STRAFTER(STR(?alt), \"http://sparqling-genomics/Sequence/\") AS ?alt
  STRAFTER(STR(?gene), \"http://sparqling-genomics/Gene/\") AS ?gene
  ?ceffect
  STRAFTER(STR(?clonality), \"http://sparqling-genomics/Clonality/\") AS ?clonality
  STRAFTER(STR(?trinucleotidecontext), \"http://sparqling-genomics/NucleotideContext/\") AS ?trinucleotidecontext
  STRAFTER(STR(?type), \"http://sparqling-genomics/VariantCallType/\") AS ?type
  STRAFTER(STR(?filter), \"http://sparqling-genomics/FilterItem/\") AS ?filter
  ?adjustedvaf
  ?adjustedcopynumber
  ?minoralleleploidy
  
  FROM <http://hmfpatients/somaticvariant>
    WHERE 
  { ?row col:sampleid ?sampleid .
    ?row col:chromosome ?chromosome .
    ?row col:position ?position .
    ?row col:ref ?ref .
    ?row col:alt ?alt .
    ?row col:type ?type .
    ?row col:gene ?gene .
    ?row col:canonicaleffect ?ceffect .
    ?row col:clonality ?clonality .
    ?row col:adjustedvaf ?adjustedvaf .
    ?row col:adjustedcopynumber ?adjustedcopynumber .
    ?row col:minoralleleploidy ?minoralleleploidy .
    ?row col:trinucleotidecontext ?trinucleotidecontext .
    ?row col:filter ?filter .
    
    FILTER (?trinucleotidecontext = trinucleotidecontext_pf:CTT AND ?ref = seq_pf:T AND ?alt = seq_pf:G ||
              ?trinucleotidecontext = trinucleotidecontext_pf:AAG AND ?ref = seq_pf:A AND ?alt = seq_pf:C)
    FILTER (?filter = filter_pf:PASS)
    FILTER (?type = type_pf:SNP)
  
  }
  ")

#connect()
query_out <- SPARQL(url = endpoint, curl_args = auth_options, query = query)$results


# split <- unlist(strsplit(sampleid, "/"))[5]
# sample_name <- substr(split, 1, nchar(split)-1)
sample_name <- sampleid



#library(devtools)
#install_github("im3sanger/dndscv")
library("dndscv")
