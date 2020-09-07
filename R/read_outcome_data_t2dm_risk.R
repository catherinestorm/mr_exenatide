

### load outcome data
  
out <- read_outcome_data(snps = exp$SNP,
                         filename = "~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/xue2018_t2dm_risk_with_eqtl_rsid.txt",
                         sep = "\t",
                         snp_col = str_c("rsid_", EXPOSURE_DATA),
                         beta_col = "b",
                         se_col = "se",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         eaf_col = "frq_A1",
                         pval_col = "P",
                         samplesize_col = "N")


out$outcome <- "T2DM risk"
  