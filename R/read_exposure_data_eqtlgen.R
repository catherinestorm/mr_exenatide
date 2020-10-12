# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")
library("readr")
source("~/Documents/UCL_PhD/PhD_Project/mr_resources/papers_resources/methods/mr_fclr/FCLR.R")


# load exposure data

exp0 <- read_exposure_data(
  filename = "data/eqtlgen_mr_exposure_dat_glp1r_pathway.txt",
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "AssessedAllele",
  other_allele_col = "OtherAllele",
  pval_col = "Pvalue",
  phenotype_col = "Gene",
  samplesize_col = "NrSamples", #double check if this is it
  min_pval = 1e-400
)

exp0 <- subset(exp0, exp0$pval.exposure < 5e-8)

exp_to_keep <- read_csv("data/genes_for_analysis.csv")
exp <- subset(exp0, (exp0$exposure %in% exp_to_keep$ensembl))

