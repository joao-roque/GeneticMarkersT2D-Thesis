
# this function returns the respective gene for each SNP
chr_to_gene <- function(list_of_variants){
  
  print('From R: OKAY')
  #source("https://bioconductor.org/biocLite.R")
  #biocLite(c("XML","biomaRt"))
  #install.packages("xml2")
  library("biomaRt")
  library("xml2")
  
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  
  # example
  #chr.region <- c("11:56468368:56468368","16:77462979:77500915")
  
  # parsing
  print('Parsing...')
  list_of_variants <- gsub("chr", "", list_of_variants)
  list_of_variants <- parse_variants_name(list_of_variants)
  
  # getting gene
  genes_list = c()
  variant_name = c()
  length_of_variants = length(list_of_variants)
  print(length_of_variants)
  print('From R: Getting genes...')
  for (i in 1:length_of_variants){
    
    tryCatch({
      filterlist <- list(list_of_variants[i],"protein_coding")
      results <- getBM(attributes = c("hgnc_symbol","entrezgene", "chromosome_name", "start_position", "end_position"), 
                    filters=c("chromosomal_region","biotype"), 
                    values = filterlist, 
                    mart = ensembl)
      
      if (length(results$hgnc_symbol) == 0){
        genes_list <- c(genes_list, "None")
        variant_name <- c(variant_name, list_of_variants[i])
      }
      else{
        # In case it returns more than 1 gene...
        if (length(results$hgnc_symbol) > 1){
            genes_list <- c(genes_list, results$hgnc_symbol[1])
            variant_name <- c(variant_name, list_of_variants[i])
        }
        else{
            genes_list <- c(genes_list, results$hgnc_symbol)
            variant_name <- c(variant_name, list_of_variants[i])
            
        }
      }
    },
    error = function(e){
      print(e)
      genes_list <- c(genes_list, "None")
      variant_name <- c(variant_name, list_of_variants[i])
      
    })

    cat("\r", i/length_of_variants)
  }
  
  newList <- list("genes" = genes_list, "variants" = variant_name)
  #genes_list
  return(newList)
  
}

parse_variants_name <- function(list_of_variants){
  
  parsed_variants <- strsplit(list_of_variants, ":")
  
  list_of_variants = c()
  for (i in 1:length(parsed_variants)){
    
    # example:
    # 11:56468368:56468368
    list_of_variants <- c(list_of_variants, 
                          paste(parsed_variants[[i]][1],
                          ":",
                          parsed_variants[[i]][2],
                          ":",
                          parsed_variants[[i]][2],
                          sep = '')
                          )
    
    
  }
  
  return(list_of_variants)
  
}
