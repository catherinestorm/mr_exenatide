# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")
library("readr")
source("~/Documents/UCL_PhD/PhD_Project/mr_resources/papers_resources/methods/mr_fclr/FCLR.R")

# load exposure data

exp0 <- read_exposure_data(
  filename = "data/psychencode_mr_exposure_dat_glp1r_pathway.txt",
  sep = "\t",
  snp_col = "Rsid",
  beta_col = "regression_slope",
  se_col = "se",
  effect_allele_col = "ALT",
  other_allele_col = "REF",
  pval_col = "nominal_pval",
  phenotype_col = "gene_id1",
  min_pval = 1e-400
)

# add sample size

exp0$samplesize.exposure <- 1387

exp0 <- subset(exp0, exp0$pval.exposure < 5e-8)

exp_to_keep <- read_csv("data/genes_for_analysis.csv")
exp <- subset(exp0, (exp0$exposure %in% exp_to_keep$ensembl))