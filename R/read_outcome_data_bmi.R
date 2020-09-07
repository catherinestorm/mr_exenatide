

### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                         filename = "data/BMI_Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",
                         sep = "\t",
                         snp_col = "SNP",
                         beta_col = "BETA",
                         se_col = "SE",
                         effect_allele_col = "Tested_Allele",
                         other_allele_col = "Other_Allele",
                         eaf_col = "Freq_Tested_Allele_in_HRS",
                         pval_col = "P",
                         samplesize_col = "N")
out$outcome <- "BMI"
