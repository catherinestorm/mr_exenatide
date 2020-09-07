library(readr)

setwd("~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020")


# read in data
eqtlgen <- read_tsv("~/Documents/UCL_PhD/PhD_Project/eqtl/eqtlgen/eqtlgen_mr_exposure_dat.txt")

psychencode <- read_tsv("~/Documents/UCL_PhD/PhD_Project/eqtl/psychencode/psychencode_mr_exposure_dat.txt")

psychencode$gene_id1 <- gsub('\\..*', '', psychencode$gene_id)

exp_to_keep <- read_csv("genes_for_analysis.csv")


# subset to keep genes in glp1r pathway
eqtlgen_glp1r_pathway <- subset(eqtlgen, (eqtlgen$Gene %in% exp_to_keep$ensembl))

psychencode_glp1r_pathway <- subset(psychencode, (psychencode$gene_id1 %in% exp_to_keep$ensembl))


# write out new files

write.table(eqtlgen_glp1r_pathway, "data/eqtlgen_mr_exposure_dat_glp1r_pathway.txt", sep = "\t", row.names = F)
write.table(psychencode_glp1r_pathway, "data/psychencode_mr_exposure_dat_glp1r_pathway.txt", sep = "\t", row.names = F)
