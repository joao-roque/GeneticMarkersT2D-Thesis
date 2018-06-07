library("biomaRt")

this.dir = dirname(parent.frame(2)$ofile)
setwd(this.dir)

variantsHomer = read.csv(file="../../data/knownvariants/homer.csv", header=TRUE, sep=",")
variantsJasper = read.csv(file="../../data/knownvariants/jasper.csv", header=TRUE, sep=",")

snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

snp_ids_homer = variantsHomer["rsID"]
snp_ids_jasper = variantsJasper["rs"]
snp_ids_jasper = unique(snp_ids_jasper)

# don't use repeated variants 
snp_ids_homer = setdiff(snp_ids_homer, snp_ids_jasper)
snp_ids_jasper = setdiff(snp_ids_jasper, snp_ids_homer)
colnames(snp_ids_jasper)[1] = "rsID"

snp_ids = merge(snp_ids_homer, snp_ids_jasper, by = "rsID")
snp_ids = unique(snp_ids)

snp_attributes = c("refsnp_id", "chr_name", "chrom_start")
snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", 
                      values=snp_ids, mart=snp_mart)

# save
write.csv(file = "../../data/knownvariants/all.csv", x = snp_locations)


