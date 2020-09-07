

### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                                     filename = str_c("~/Documents/UCL_PhD/PhD_Project/pdgwas_data/progression/",OUTCOME,"_with_eqtl_rsid.txt"),
                                     sep = ",",
                                     snp_col = str_c("rsid_", EXPOSURE_DATA),
                                     beta_col = "BETA",
                                     se_col = "SE",
                                     effect_allele_col = "effect_allele",
                                     other_allele_col = "other_allele",
                                     pval_col = "P",
                                     samplesize_col = "N",
                                     phenotype_col = "outcome"
                                     )
