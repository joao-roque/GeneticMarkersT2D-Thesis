chr_to_rs <- function(list_of_variants){
  
  print('From R: OKAY')
  #source("https://bioconductor.org/biocLite.R")
  #biocLite(c("XML","biomaRt"))
  #install.packages("xml2")
  library("biomaRt")
  library("xml2")
  
  #ensembl <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_gene_ensembl")
  ensembl <- useMart("ENSEMBL_MART_SNP", "hsapiens_snp")
  
  # example
  #chr.region <- c("11:56468368:56468368","16:77462979:77500915")
  
  # parsing
  print('Parsing...')
  list_of_variants <- gsub("chr", "", list_of_variants)
  list_of_variants <- parse_variants_name(list_of_variants)
  
  #list_of_variants <- c("11:56468368:56468368","16:77462979:77500915")
  
  # getting gene
  rs_list = c()
  variant_name = c()
  length_of_variants = length(list_of_variants)
  #print(length_of_variants)
  print('From R: Getting genes...')
  
  #print(list_of_variants)
  
  results <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"), 
                   filters = "chromosomal_region", 
                   values = list_of_variants, 
                   mart = ensembl)
  
  #for (i in 1:length_of_variants){
    # tryCatch({
    #   
    #   #variant = strsplit(list_of_variants[i],":")
    #   #filterlist <- list(variant[[1]][1], variant[[1]][2], variant[[1]][3])#, "protein_coding")
    #   
    #   results <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start", "chrom_end"), 
    #                    filters = "chromosomal_region", 
    #                    values = list_of_variants, 
    #                    mart = ensembl)
    #   
    #   if (length(results$refsnp_id) == 0){
    #     rs_list <- c(rs_list, "None")
    #     variant_name <- c(variant_name, list_of_variants[i])
    #   }
    #   else{
    #     # In case it returns more than 1 gene...
    #     if (length(results$refsnp_id) > 1){
    #       rs_list <- c(rs_list, results$refsnp_id[1])
    #       variant_name <- c(variant_name, list_of_variants[i])
    #     }
    #     else{
    #       rs_list <- c(rs_list, results$refsnp_id)
    #       variant_name <- c(variant_name, list_of_variants[i])
    #       
    #     }
    #   }
    # },
    # error = function(e){
    #   print(e)
    #   rs_list <- c(rs_list, "None")
    #   variant_name <- c(variant_name, list_of_variants[i])
    #   
    # })
      
      
      #cat("\r", i/length_of_variants)
    #}
  
  
  
  
  newList <- list("rsID" = results)
  write.csv(file = "../../data/knownvariants/rsID.csv", x =newList)
  #genes_list
  return(newList)
  
}

parse_variants_name <- function(list_of_variants){
  
  parsed_variants <- strsplit(list_of_variants, ":")

  list_of_variants = c()
  for (i in 1:length(parsed_variants)){
    
    # example:
    # 11:56468368:56468368
    # variant_plus_one <-  toString(as.numeric(parsed_variants[[i]][2]) + 1)
    variant <- toString(parsed_variants[[i]][2])
    
    list_of_variants <- c(list_of_variants, 
                          paste(parsed_variants[[i]][1],
                                ":",
                                parsed_variants[[i]][2],
                                ":",
                                variant,
                                sep = '')
    )
    
    
  }
  
  return(list_of_variants)
  
}
