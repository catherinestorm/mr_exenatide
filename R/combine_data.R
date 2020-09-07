setwd("~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020/results_5e8")
library(readr)
library(dplyr)


exp_to_keep <- read_csv("~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020/genes_for_analysis.csv")





# CLASSIC MR METHODS r2 0.2
classic_files_0.2 <- list.files(pattern = "results_liberal_r2_0.2")
classic_files_0.2 <- classic_files_0.2[!grepl("egger",classic_files_0.2)]
classic_files_0.2 <- classic_files_0.2[!grepl("steiger",classic_files_0.2)]

all_res_0.2 <- data.frame()

for (i in 1:length(classic_files_0.2)) {
  temp <- read_tsv(classic_files_0.2[i])
  temp <- subset(temp, temp$prot_gene %in% exp_to_keep$prot_gene)
  temp$tissue <- classic_files_0.2[i]
  temp$tissue <- gsub("results_liberal_r2_0.2_", "", temp$tissue)
  temp$tissue <- gsub("_.*", "", temp$tissue)
  all_res_0.2 <- distinct(rbind(all_res_0.2, temp))
}

all_res_0.2 <- subset(all_res_0.2, !(all_res_0.2$outcome == "T2DM risk" & all_res_0.2$tissue == "psychencode"))

data_for_qval <- which(all_res_0.2$method == "IVW" | all_res_0.2$method == "Inverse variance weighted" | all_res_0.2$method == "Wald ratio")
all_res_0.2$fdr_new <- NA
all_res_0.2$fdr_new[data_for_qval] <- p.adjust(all_res_0.2$p[data_for_qval], method = "fdr")

all_res_0.2_sign <- subset(all_res_0.2, all_res_0.2$fdr_qval < 0.05)
all_res_0.2_sign <- subset(all_res_0.2, all_res_0.2$fdr_new < 0.05)

bonferroni <- 0.05/(length(unique(all_res_0.2$prot_gene))*length(unique(all_res_0.2$outcome)))
bonferroni <- 0.05/18

all_res_0.2_sign <- subset(all_res_0.2, all_res_0.2$p < bonferroni)

write.table(all_res_0.2, "full_res_liberal_r2_0.2.txt", row.names = F, sep = "\t")



# STEIGER
classic_files_0.2_steiger <- list.files(pattern = "steiger_results_liberal_r2_0.2")

all_res_0.2_steiger <- data.frame()

for (i in 1:length(classic_files_0.2_steiger)) {
  temp <- read_tsv(classic_files_0.2_steiger[i])
  temp$tissue <- classic_files_0.2_steiger[i]
  temp$tissue <- gsub("steiger_results_liberal_r2_0.2_", "", temp$tissue)
  temp$tissue <- gsub("_.*", "", temp$tissue)
  all_res_0.2_steiger <- distinct(rbind(all_res_0.2_steiger, temp))
}

all_res_0.2_steiger <- subset(all_res_0.2_steiger, !(all_res_0.2_steiger$outcome == "T2DM risk" & all_res_0.2_steiger$tissue == "psychencode"))

all_res_0.2_steiger <- left_join(exp_to_keep,all_res_0.2_steiger, by = c("ensembl" = "exposure"))


write.table(all_res_0.2_steiger, "full_steiger_full_results_liberal_r2_0.2.txt", row.names = F, sep = "\t")






# CLASSIC MR METHODS r2 0.4
classic_files_0.4 <- list.files(pattern = "results_liberal_r2_0.4")
classic_files_0.4 <- classic_files_0.4[!grepl("egger",classic_files_0.4)]

all_res_0.4 <- data.frame()

for (i in 1:length(classic_files_0.4)) {
  temp <- read_tsv(classic_files_0.4[i])
  temp <- subset(temp, temp$prot_gene %in% exp_to_keep$prot_gene)
  temp$tissue <- classic_files_0.4[i]
  temp$tissue <- gsub("results_liberal_r2_0.4_", "", temp$tissue)
  temp$tissue <- gsub("_.*", "", temp$tissue)
  all_res_0.4 <- distinct(rbind(all_res_0.4, temp))
}


all_res_0.4_sign <- subset(all_res_0.4, all_res_0.4$fdr_qval < 0.05)


# CLASSIC MR METHODS r2 0.001
classic_files_strict <- list.files(pattern = "results_liberal_r2_0.001")

all_res_strict <- data.frame()

for (i in 1:length(classic_files_strict)) {
  temp <- read_tsv(classic_files_strict[i])
  temp <- subset(temp, temp$prot_gene %in% exp_to_keep$prot_gene)
  all_res_strict <- distinct(rbind(all_res_strict, temp))
}

names(all_res_strict)[names(all_res_strict) == "exposure_data"] <- "tissue"

all_res_strict_sign <- subset(all_res_strict, all_res_strict$fdr_qval < 0.05)


all_clump <- distinct(rbind(all_res_0.2, all_res_strict[,names(all_res_0.2)]))

all_clump <- distinct(rbind(all_clump, all_res_0.4[,names(all_clump)]))

write.table(all_clump, "full_results_mr_classic.txt", row.names = F, sep = "\t")




# bonferroni

bonferroni_threshold <- 0.05/length(unique(all_clump$prot_gene))


all_clump$prot_gene_tissue <- paste(all_clump$prot_gene, all_clump$tissue, sep = "_")

bonferroni_threshold <- 0.05/length(unique(all_clump$prot_gene_tissue))



all_clump$prot_gene_tissue_outcome <- paste(all_clump$prot_gene, all_clump$tissue, all_clump$outcome, sep = "_")

bonferroni_threshold <- 0.05/length(unique(all_clump$prot_gene_tissue_outcome))


all_clump_bonf_sign <- subset(all_clump, all_clump$p < bonferroni_threshold)

all_clump_fdr_sign <- subset(all_clump, all_clump$fdr_qval < 0.05)




# glp1r, dpp4, tlr4

all_clump_selected <- subset(all_clump, all_clump$prot_gene %in% c("TLR4","GLP1R","DPP4"))

all_clump_selected <- subset(all_clump_selected, all_clump_selected$method %in% c("Inverse variance weighted", "IVW", "Wald ratio") & all_clump_selected$clump_thres == 0.2)
all_clump_selected$fdr_new <- p.adjust(all_clump_selected$p, method = "fdr")

bonferroni_threshold_select <- 0.05/length(unique(all_clump_selected$prot_gene_tissue_outcome))

all_clump_selected_bonf_sign <- subset(all_clump_selected, all_clump_selected$p < bonferroni_threshold_select)




# QCMETHODS r2 0.2
qc_files <- list.files(pattern = "results_liberal_r2_0.2")

all_qc <- data.frame()

for (i in 1:length(qc_files)) {
  temp <- read_tsv(qc_files[i])
  temp <- subset(temp, temp$prot_gene %in% exp_to_keep$prot_gene)
  temp$tissue <- qc_files[i]
  temp$tissue <- gsub("results_liberal_r2_0.2_egger_cochransq_", "", temp$tissue)
  temp$tissue <- gsub("_.*", "", temp$tissue)
  all_qc <- distinct(rbind(all_qc, temp))
}

all_qc_select <- subset(all_qc, all_qc$prot_gene %in% c("GLP1R","TLR4","DPP4"))

write.table(all_qc, "full_results_mr_classic_qc.txt", row.names = F, sep = "\t")







# PCA
mr_pca_files <- list.files(pattern = "results_pca")

mr_pca <-  data.frame()

for (i in 1:length(mr_pca_files)) {
  temp <- read_tsv(mr_pca_files[i])
  temp <- subset(temp, temp$prot_gene %in% exp_to_keep$prot_gene)
  temp$fdr_qval <- NA
  temp$fdr_qval <- p.adjust(temp$p, method = "fdr")
  
  mr_pca <- distinct(rbind(mr_pca, temp))
}


mr_pca_sign <- mr_pca[mr_pca$fdr_qval < 0.05,]

write.table(mr_pca, "full_results_mr_pca.txt", row.names = F, sep = "\t")






# POWER CALCULATION

expit <- function(x) { return(exp(x)/(1+exp(x))) }
rsq = 0.005547943 # squared correlation

b1 = -0.061650440 # causal effect (log odds ratio per SD
#b1 = log(1.2) # or log of OR per SD)

sig = 0.05 # significance level (alpha)
pow = 0.8 # power level (1-beta)
ratio = 1 # ratio of cases:controls = 1:ratio

cat("Sample size required for ", pow*100, "% power: ",
    (qnorm(1-sig/2)+qnorm(pow))^2/b1^2/rsq/(ratio/(1+ratio))/(1/(1+ratio)))

n = 485648 # Sample size
cat("Power of analysis with ", n, "participants: ",
    pnorm(sqrt(n*rsq*(ratio/(1+ratio))*(1/(1+ratio)))*b1-qnorm(1-sig/2)))
