


library("dplyr")
library("readr")
library("stringr")

## read in eqtl data

eqtlgen <- read.table("data/eqtlgen_mr_exposure_dat_glp1r_pathway.txt", sep = "\t",colClasses = "character", header = T)

eqtlgen <- distinct(eqtlgen[,c("SNP", "SNPChr", "SNPPos")])

names(eqtlgen) <- c("rsid_eqtlgen", "chr", "position")

eqtlgen$chr_pos <- str_c("chr",eqtlgen$chr, ":",eqtlgen$position)


psychencode <- read.table("data/psychencode_mr_exposure_dat_glp1r_pathway.txt", sep = "\t",colClasses = "character", header = T)

psychencode <- distinct(psychencode[,c("Rsid", "chr", "position")])

names(psychencode) <- c("rsid_psychencode", "chr", "position")

psychencode$chr_pos <- str_c(psychencode$chr, ":",psychencode$position)





## read in pd risk meta5 gwas data
nalls2019 <- read_tsv(gzfile("data/meta_no23.tbl.gz"))

nalls2019 <- left_join(nalls2019, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = c("MarkerName" = "chr_pos"))
nalls2019 <- left_join(nalls2019, psychencode[,c("chr_pos", "rsid_psychencode")], by = c("MarkerName" = "chr_pos"))


write.table(nalls2019, "data/nalls2019_with_eqtl_rsid.txt", sep = "\t", row.names = F, quote = F)



## read in pd aao gwas data
aao <- read_tsv(unzip("data/Blauwendraat_IPDGC_only_AAO_GWAS_sumstats_april_2018.txt.zip"))

aao <- left_join(aao, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = c("MarkerName" = "chr_pos"))
aao <- left_join(aao, psychencode[,c("chr_pos", "rsid_psychencode")], by = c("MarkerName" = "chr_pos"))


write.table(aao, "data/blauwendraat2019_aao_with_eqtl_rsid.txt", sep = "\t", row.names = F, quote = F)





## read in pd progression data
alleles <- read.table(gzfile("data/progression/reference.txt.gz"), sep = ",", stringsAsFactors = F, header = T)
names(alleles) <- c("chrpos", "rsid", "chr", "start", "other_allele", "effect_allele", "maf", "func", "near_gene")


file_list_res_a <- list.files("data/progression", pattern = "cont_", full.names = T)
file_list_res_b <- list.files("data/progression", pattern = "surv_", full.names = T)
file_list_res <- c(file_list_res_a, file_list_res_b)


for (i in 1:length(file_list_res)) {
  raw_data <- read.table(gzfile(file_list_res[i]), sep = "\t", stringsAsFactors = F, header = T)
  
  raw_data$outcome <- file_list_res[i]
  raw_data$outcome <- gsub(".txt.gz", "", raw_data$outcome)
  raw_data$outcome <- gsub("data/progression/", "", raw_data$outcome)
  
  raw_data$SNP1 <- str_c("chr", raw_data$SNP)
  
  raw_data_with_alleles <- left_join(raw_data, alleles[, c("chrpos", "rsid", "other_allele", "effect_allele", "maf")], by = c("SNP" = "chrpos"))
  
  raw_data_with_alleles_1 <- left_join(raw_data_with_alleles, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = c("SNP1" = "chr_pos"))
  raw_data_with_alleles_1 <- left_join(raw_data_with_alleles_1, psychencode[,c("chr_pos", "rsid_psychencode")], by = c("SNP1" = "chr_pos"))
  
  
  write.table(raw_data_with_alleles_1, str_c("data/progression/", paste(raw_data_with_alleles$outcome[1]), "_with_eqtl_rsid.txt"), sep = ",", row.names = F, quote = T)
}




## read in t2dm risk gwas data
t2dm <- read_csv("data/Xue_et_al_T2D_META_Nat_Commun_2018.txt")
t2dm$chr_pos <- str_c("chr", t2dm$CHR, ":", t2dm$BP)

t2dm <- left_join(t2dm, eqtlgen[,c("chr_pos", "rsid_eqtlgen")], by = "chr_pos")
t2dm <- left_join(t2dm, psychencode[,c("chr_pos", "rsid_psychencode")], by = "chr_pos")

write.table(t2dm, "data/xue2018_t2dm_risk_with_eqtl_rsid.txt", sep = "\t", row.names = F, quote = F)

