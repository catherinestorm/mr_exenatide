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
all_res_0.2 <- subset(all_res_0.2, !(all_res_0.2$outcome == "BMI" & all_res_0.2$tissue == "psychencode"))

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
all_res_0.2_steiger <- subset(all_res_0.2_steiger, !(all_res_0.2_steiger$outcome == "BMI" & all_res_0.2_steiger$tissue == "psychencode"))

all_res_0.2_steiger <- left_join(exp_to_keep,all_res_0.2_steiger, by = c("ensembl" = "exposure"))


write.table(all_res_0.2_steiger, "full_steiger_full_results_liberal_r2_0.2.txt", row.names = F, sep = "\t")





# QCMETHODS r2 0.2
qc_files <- list.files(pattern = "results_liberal_r2_0.2_egger_cochransq")

all_qc <- data.frame()

for (i in 1:length(qc_files)) {
  temp <- read_tsv(qc_files[i])
  temp <- subset(temp, temp$prot_gene %in% exp_to_keep$prot_gene)
  temp$tissue <- qc_files[i]
  temp$tissue <- gsub("results_liberal_r2_0.2_egger_cochransq_", "", temp$tissue)
  temp$tissue <- gsub("_.*", "", temp$tissue)
  all_qc <- distinct(rbind(all_qc, temp))
}

all_qc <- subset(all_qc, !(all_qc$outcome == "T2DM risk" & all_qc$tissue == "psychencode"))
all_qc <- subset(all_qc, !(all_qc$outcome == "BMI" & all_qc$tissue == "psychencode"))


all_qc_select <- subset(all_qc, all_qc$prot_gene %in% c("GLP1R","TLR4","DPP4"))

write.table(all_qc, "full_results_mr_classic_qc.txt", row.names = F, sep = "\t")




