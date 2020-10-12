library(readr)


# EQTLGEN

# load eqtlgen data
eqtlgen0 <- read_tsv(gzfile("data/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"), col_types=cols(.default="c"))

eqtlgen0$new_gene_id <- eqtlgen0$Gene


# alt allele == effect allele
alleles <- read_tsv(gzfile("data/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"), col_types=cols(.default="c"))


# AlleleB_all ==  Allele frequency of Allele B
full <- left_join(eqtlgen0, alleles[, c("SNP", "AlleleB", "AlleleB_all")], by = "SNP")


#switch allele frequencies where appropriate
mismatch <- which(full$AssessedAllele != full$AlleleB)
mismatch2 <- which(full$OtherAllele == full$AlleleB)
sum(mismatch - mismatch2) # should be 0

full$eaf <- full$AlleleB_all
full$eaf[mismatch] <- 1 - as.numeric(full$eaf[mismatch])



# calculate beta and standard error

full$beta <- as.numeric(full$Zscore) / sqrt(2 * as.numeric(full$eaf) * 
                                              (1 - as.numeric(full$eaf)) * 
                                              (as.numeric(full$NrSamples) + as.numeric(full$Zscore)^2))

full$se <- 1 / sqrt(2 * as.numeric(full$eaf) * 
                      (1 - as.numeric(full$eaf)) * 
                      (as.numeric(full$NrSamples) + as.numeric(full$Zscore)^2))

eqtlgen <- full 








# PSYCHENCODE

psychencode0 <- read_tsv("psychencode/DER-08a_hg19_eQTL.significant.txt", col_types=cols(.default="c"))

# alt allele == effect allele
alleles <- read_tsv("psychencode/SNP_Information_Table_with_Alleles.txt", col_types=cols(.default="c"))


# add allele data
full <- left_join(psychencode0, alleles, by = c("SNP_id" = "PEC_id"))


# calculate standard error from beta and pvalue

full$se=abs(as.numeric(full$regression_slope)/qnorm(as.numeric(full$nominal_pval)/2))

all(full$SNP_start == full$SNP_end) # TRUE

psychencode <- full

psychencode$gene_id1 <- gsub('\\..*', '', psychencode$gene_id)




# SUBSET TO KEEP GENES IN GLP1R PATHWAY ONLY

exp_to_keep <- read_csv("genes_for_analysis.csv")

eqtlgen_glp1r_pathway <- subset(eqtlgen, (eqtlgen$Gene %in% exp_to_keep$ensembl))

psychencode_glp1r_pathway <- subset(psychencode, (psychencode$gene_id1 %in% exp_to_keep$ensembl))


# write out new files

write.table(eqtlgen_glp1r_pathway, "data/eqtlgen_mr_exposure_dat_glp1r_pathway.txt", sep = "\t", row.names = F)
write.table(psychencode_glp1r_pathway, "data/psychencode_mr_exposure_dat_glp1r_pathway.txt", sep = "\t", row.names = F)
