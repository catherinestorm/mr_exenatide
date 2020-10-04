setwd("~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020/results")

# prep
library(dplyr)
library(readr)
library(forestplot)
library(stringr)
library(tidyverse)

# this function will format it's input to scientific notation if the number is < 0.001 and round it to 3 significant figures if it's >= 0.001
# useful for p values
format_numbers <- function(number) {
  if (is.numeric(number) == FALSE & is.integer(number) == FALSE) {
    result <- NA
  #} else if (number < 0.0001) {
    #result <- formatC(number, format = "e", digits = 2)
  } else { #if (number >= 0.0001) {
    result <- format(round(number, 3), nsmall = 3)
  }
  return(result)
}




# make forest plot T2DM--------

data0 <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_t2dm_risk_all.txt")) # read in data
data0$tissue <- "eqtlgen"

data <- subset(data0, data0$prot_gene %in% c("GLP1R","TLR4","DPP4"))

# calculate odds ratio and sort by gene names
data <- TwoSampleMR::generate_odds_ratios(data)
#data <- data[order(data$prot_gene),]
data <- data[order(match(data$prot_gene, c("GLP1R","DPP4", "TLR4"))),]

# subset to keep IVW and Wald ratio only
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]

# tissue names if diff. tissues
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"

# create a column with ODDS RATIO (CONFIDENCE INTERVAL) and format the pvalues
data$or_ci <- str_c(format_numbers(data$or), " (", format_numbers(data$or_lci95), ", ",format_numbers(data$or_uci95), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

# DATA FOR THE CREATING THE FOREST PLOT
forest_data <- data[,c("or", "or_lci95", "or_uci95", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data) # size we have a header row we want a blank ros at the top

# DATA TO DISPLAY AS NUMBETS IM THE THE FOREST PLOT
table_text <- data[, c("prot_gene","nsnp", "or_ci", "roundp")]
table_text <- rbind(c("Gene", "No. SNPs", "OR (95% CI)", "FDR-adjusted P"), table_text)


pdf("plots/forest_t2dm_risk.pdf", width = 23,height = 10, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.1, # size of the box
           align = "l", # align text to left 
           mean = forest_data$or,
           lower = forest_data$or_lci95,
           upper = forest_data$or_uci95,
           xlog=TRUE, # log scale if odds ratios
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")),  # all "Blood" rows are now "summary" rows; this is a trick to allow data from diff tissues to be colour coded
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours for blood vs brain tissue
           xlab="T2DM odds ratio",
           hrzl_lines = list("2" = gpar(lty = 2)), # a horizontal line at the top under the header
           txt_gp = fpTxtGp(summary = lapply(c(3,1,1,1), 
                                             function(val)  gpar(fontface = val)),
                            label = gpar(cex=3),
                            ticks = gpar(cex=3),
                            xlab  = gpar(fontface="bold", cex = 3)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0.1, ...)} # format the "summary" rows to suit the main rows
)
dev.off()








# make forest plot BMI----------

data0 <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_bmi_all.txt"))
data0$tissue <- "eqtlgen"

data <- subset(data0, data0$prot_gene %in% c("GLP1R","TLR4","DPP4"))

# calculate confidence intervals and sort by gene names
data <- TwoSampleMR::generate_odds_ratios(data)

#data <- data[order(data$prot_gene),]
data <- data[order(match(data$prot_gene, c("GLP1R","DPP4", "TLR4"))),]

# subset to keep IVW and Wald ratio only
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]

# tissue names if diff. tissues
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"

# create a column with BETA (CONFIDENCE INTERVAL) and format the pvalues
data$beta_ci <- str_c(format_numbers(data$beta), " (", format_numbers(data$lo_ci), ", ",format_numbers(data$up_ci), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)


# DATA FOR THE CREATING THE FOREST PLOT
forest_data <- data[,c("beta", "lo_ci", "up_ci", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

# DATA TO DISPLAY AS NUMBETS IM THE THE FOREST PLOT
table_text <- data[, c("prot_gene","nsnp", "beta_ci", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text <- rbind(c("Gene", "No. SNPs", "Beta (95% CI)", "FDR-adjusted P"), table_text)


pdf("plots/forest_bmi.pdf", width = 23,height = 10, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.1, # size of the box
           align = "l", # align text to left 
           mean = forest_data$beta,
           lower = forest_data$lo_ci,
           upper = forest_data$up_ci,
           xlog=F, # no log scale bc beta
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows; this is a trick to allow data from diff tissues to be colour coded
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours for blood vs brain tissue
           xlab="SD-change in BMI",
           hrzl_lines = list("2" = gpar(lty = 2)), # a horizontal line at the top under the header
           txt_gp = fpTxtGp(summary = (lapply(c(3,1,1,1), 
                                              function(val)  gpar(fontface = val, cex = 3))),
                            label = gpar(cex=3),
                            ticks = gpar(cex=3),
                            xlab  = gpar(fontface="bold", cex = 3)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0.1, ...)} # format the "summary" rows to suit the main rows
)
dev.off()







# make forest plot PD--------

data01 <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_pd_risk_all.txt"))
data01$tissue <- "eqtlgen"
data02 <- data.frame(read_tsv("results_liberal_r2_0.2_psychencode_pd_risk_all.txt"))
data02$tissue <- "psychencode"

data0 <- distinct(rbind(data01, data02))

data <- subset(data0, data0$prot_gene %in% c("GLP1R","TLR4","DPP4"))



data <- TwoSampleMR::generate_odds_ratios(data)
#data <- data[order(data$prot_gene),]
data <- data[order(match(data$prot_gene, c("GLP1R","DPP4", "TLR4"))),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
data$or_ci <- str_c(format_numbers(data$or), " (", format_numbers(data$or_lci95), ", ",format_numbers(data$or_uci95), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

forest_data <- data[,c("or", "or_lci95", "or_uci95", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

table_text <- data[, c("prot_gene","nsnp", "or_ci", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text <- rbind(c("Gene", "No. SNPs", "OR (95% CI)", "FDR-adjusted P"), table_text)


pdf("plots/forest_pd_risk.pdf", width = 23,height = 13, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.1,
           align = "l",
           mean = forest_data$or,
           lower = forest_data$or_lci95,
           upper = forest_data$or_uci95,
           xlog=TRUE, # log scale
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours
           xlab="Parkinson's odds ratio",
           hrzl_lines = list("2" = gpar(lty = 2)),
           txt_gp = fpTxtGp(summary = (lapply(c(3,1,1,1), 
                                             function(val)  gpar(fontface = val, cex = 3))),
                            label = (lapply(c(3,1,1,1), 
                                            function(val)  gpar(fontface = val, cex = 3))),
                            ticks = gpar(cex=3, fontface = "plain"),
                            xlab  = gpar(fontface="bold", cex = 3)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0.1, ...)}
)
dev.off()








# make forest plot PD AAO----------

data01 <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_pd_aao_all.txt"))
data01$tissue <- "eqtlgen"
data02 <- data.frame(read_tsv("results_liberal_r2_0.2_psychencode_pd_aao_all.txt"))
data02$tissue <- "psychencode"

data0 <- distinct(rbind(data01, data02))

data <- subset(data0, data0$prot_gene %in% c("GLP1R","TLR4","DPP4"))


data <- TwoSampleMR::generate_odds_ratios(data)
#data <- data[order(data$prot_gene),]
data <- data[order(match(data$prot_gene, c("GLP1R","DPP4", "TLR4"))),]

data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"
data$beta_ci <- str_c(format_numbers(data$beta), " (", format_numbers(data$lo_ci), ", ",format_numbers(data$up_ci), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

forest_data <- data[,c("beta", "lo_ci", "up_ci", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

table_text <- data[, c("prot_gene","nsnp", "beta_ci", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text <- rbind(c("Gene", "No. SNPs", "Beta (95% CI)", "FDR-adjusted P"), table_text)


pdf("plots/forest_pd_aao.pdf", width = 23,height = 13, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.1,
           align = "l",
           mean = forest_data$beta,
           lower = forest_data$lo_ci,
           upper = forest_data$up_ci,
           xlog=F, # log scale
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours
           xlab="SD-change in Parkinson's age at onset",
           hrzl_lines = list("2" = gpar(lty = 2)),
           txt_gp = fpTxtGp(summary = (lapply(c(3,1,1,1), 
                                              function(val)  gpar(fontface = val, cex = 3))),
                            label = (lapply(c(3,1,1,1), 
                                            function(val)  gpar(fontface = val, cex = 3))),
                            ticks = gpar(cex=3, fontface="plain"),
                            xlab  = gpar(fontface="bold", cex = 3)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0.1, ...)}
)
dev.off()








# MAKE A FOREST PLOT FOR ALL THE TISSUES
data_t2dm <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_t2dm_risk_all.txt"))
data_t2dm$tissue <- "eqtlgen"

data_bmi <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_bmi_all.txt"))
data_bmi$tissue <- "eqtlgen"

data_t2dm_bmi <- distinct(rbind(data_t2dm, data_bmi))

data_pd_risk_eqtlgen <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_pd_risk_all.txt"))
data_pd_risk_eqtlgen$tissue <- "eqtlgen"
data_pd_risk_psychencode <- data.frame(read_tsv("results_liberal_r2_0.2_psychencode_pd_risk_all.txt"))
data_pd_risk_psychencode$tissue <- "psychencode"

data_pd_risk <- distinct(rbind(data_pd_risk_eqtlgen, data_pd_risk_psychencode))
data_pd_risk$outcome <- "Parkinson's risk"

data_t2dm_bmi_pd_risk <- distinct(rbind(data_t2dm_bmi, data_pd_risk))

data_pd_aao_eqtlgen <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_pd_aao_all.txt"))
data_pd_aao_eqtlgen$tissue <- "eqtlgen"
data_pd_aao_psychencode <- data.frame(read_tsv("results_liberal_r2_0.2_psychencode_pd_aao_all.txt"))
data_pd_aao_psychencode$tissue <- "psychencode"

data_pd_aao <- distinct(rbind(data_pd_aao_eqtlgen, data_pd_aao_psychencode))
data_pd_aao$outcome <- "Parkinson's age at onset"

data_all_outcomes <- distinct(rbind(data_t2dm_bmi_pd_risk, data_pd_aao))

#data <- subset(data_all_outcomes, data_all_outcomes$prot_gene %in% c("GLP1R","TLR4","DPP4"))

# for t2dm, bmi, pd risk we only want GLP1R and DPP4; for AAO, we only want TlR4 in brain tissue
data <- subset(data_all_outcomes, (data_all_outcomes$prot_gene %in% c("GLP1R","DPP4") & 
                                     !(data_all_outcomes$outcome == "Parkinson's age at onset")) 
               | (data_all_outcomes$prot_gene == "TLR4" & 
                    data_all_outcomes$outcome == "Parkinson's age at onset" & data_all_outcomes$tissue == "eqtlgen"))


# calculate confidence intervals and put in desired order
data <- TwoSampleMR::generate_odds_ratios(data)
data <- data[order(match(data$prot_gene, c("GLP1R","DPP4", "TLR4"))),]

# subset to keep IVW and Wald ratio only
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]

# tissue names if diff. tissues
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"

# create a column with BETA (CONFIDENCE INTERVAL) and format the pvalues
data$beta_ci <- str_c(round(data$beta, digits = 3), " (", round(data$lo_ci, digits = 2), ", ",round(data$up_ci, digits = 2), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

# DATA FOR THE CREATING THE FOREST PLOT
forest_data <- data[,c("beta", "lo_ci", "up_ci", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

# DATA TO DISPLAY AS NUMBETS IM THE THE FOREST PLOT
table_text <- data[, c("prot_gene","outcome", "nsnp", "beta_ci", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text <- rbind(c("Gene","Outcome",  "No. SNPs", "Beta (95% CI)", "FDR-adjusted P"), table_text)


pdf("plots/forest_all.pdf", width = 30,height = 15, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.1, # size of the box
           #ci.vertices = T,
           #lty.lwd = 10,
           align = "l",
           mean = forest_data$beta,
           lower = forest_data$lo_ci,
           upper = forest_data$up_ci,
           xlog=F, # no log scale bc beta
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows; this is a trick to allow data from diff tissues to be colour coded
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours for blood vs brain tissue
           xlab="SD-change in outcome",
           clip =c(-0.5, 0.8),
           hrzl_lines = list("2" = gpar(lty = 4), # a horizontal line at the top under the header
                             "5" = gpar(lty = 2), # a horizontal line to divide the plot by gene
                             "9" = gpar(lty = 2)), # a horizontal line to divide the plot by gene
           txt_gp = fpTxtGp(summary = list(gpar(cex=2.5, fontface = "plain")), # text size etc
                            label = gpar(cex=2.5),
                            ticks = gpar(cex=2),
                            xlab  = gpar(fontface="bold", cex = 2)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0.1, ...)} # format the "summary" rows to suit the main rows
           )
dev.off()















# FOR THE POSTER MAKE A FOREST PLOT FOR ALL THE TISSUES
data_t2dm <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_t2dm_risk_all.txt"))
data_t2dm$tissue <- "eqtlgen"

data_bmi <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_bmi_all.txt"))
data_bmi$tissue <- "eqtlgen"

data_t2dm_bmi <- distinct(rbind(data_t2dm, data_bmi))

data_pd_risk_eqtlgen <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_pd_risk_all.txt"))
data_pd_risk_eqtlgen$tissue <- "eqtlgen"
data_pd_risk_psychencode <- data.frame(read_tsv("results_liberal_r2_0.2_psychencode_pd_risk_all.txt"))
data_pd_risk_psychencode$tissue <- "psychencode"

data_pd_risk <- distinct(rbind(data_pd_risk_eqtlgen, data_pd_risk_psychencode))
data_pd_risk$outcome <- "Parkinson's risk"

data_t2dm_bmi_pd_risk <- distinct(rbind(data_t2dm_bmi, data_pd_risk))

data_pd_aao_eqtlgen <- data.frame(read_tsv("results_liberal_r2_0.2_eqtlgen_pd_aao_all.txt"))
data_pd_aao_eqtlgen$tissue <- "eqtlgen"
data_pd_aao_psychencode <- data.frame(read_tsv("results_liberal_r2_0.2_psychencode_pd_aao_all.txt"))
data_pd_aao_psychencode$tissue <- "psychencode"

data_pd_aao <- distinct(rbind(data_pd_aao_eqtlgen, data_pd_aao_psychencode))
data_pd_aao$outcome <- "Parkinson's age at onset"

data_all_outcomes <- distinct(rbind(data_t2dm_bmi_pd_risk, data_pd_aao))

#data <- subset(data_all_outcomes, data_all_outcomes$prot_gene %in% c("GLP1R","TLR4","DPP4"))

# for t2dm, bmi, pd risk we only want GLP1R and DPP4; for AAO, we only want TlR4 in brain tissue
data <- subset(data_all_outcomes, (data_all_outcomes$prot_gene %in% c("GLP1R","DPP4") & 
                                     !(data_all_outcomes$outcome == "Parkinson's age at onset")) 
               | (data_all_outcomes$prot_gene == "TLR4" & 
                    data_all_outcomes$outcome == "Parkinson's age at onset" & data_all_outcomes$tissue == "eqtlgen"))


# calculate confidence intervals and put in desired order
data <- TwoSampleMR::generate_odds_ratios(data)
data <- data[order(match(data$prot_gene, c("GLP1R","DPP4", "TLR4"))),]

# subset to keep IVW and Wald ratio only
data <- data[which(data$method == "IVW" | data$method == "Inverse variance weighted" | data$method == "Wald ratio"),]

# tissue names if diff. tissues
data[which(data$tissue == "eqtlgen"), "tissue"] <- "Blood"
data[which(data$tissue == "psychencode"), "tissue"] <- "Brain"

# create a column with BETA (CONFIDENCE INTERVAL) and format the pvalues
data$beta_ci <- str_c(round(data$beta, digits = 3), " (", round(data$lo_ci, digits = 2), ", ",round(data$up_ci, digits = 2), ")")
data$roundp <- lapply(X = data$fdr_qval, FUN = format_numbers)

# DATA FOR THE CREATING THE FOREST PLOT
forest_data <- data[,c("beta", "lo_ci", "up_ci", "tissue")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data)

# DATA TO DISPLAY AS NUMBETS IM THE THE FOREST PLOT
table_text <- data[, c("prot_gene","outcome", "roundp")]
#table_text <- table_text[complete.cases(table_text),]
table_text <- rbind(c("Gene","Outcome",  "FDR-adjusted P"), table_text)


pdf("/Users/catherinestorm/Documents/UCL_PhD/PhD_Project/submissions_courses_conferences/conference_parkinsons_uk_2020/forest_all.pdf", width = 7,height = 5, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.15, # size of the box
           #ci.vertices = T,
           #lty.lwd = 10,
           align = "l",
           mean = forest_data$beta,
           lower = forest_data$lo_ci,
           upper = forest_data$up_ci,
           xlog=F, # no log scale bc beta
           is.summary= c(forest_data$tissue %in% c("Blood", "NA")), # all "Blood" rows are now "summary" rows; this is a trick to allow data from diff tissues to be colour coded
           col=fpColors(box="royalblue", line = "royalblue", summary = "darkred"), # colours for blood vs brain tissue
           xlab="SD-change in outcome",
           clip =c(-0.5, 0.8),
           hrzl_lines = list("2" = gpar(lty = 4), # a horizontal line at the top under the header
                             "5" = gpar(lty = 2), # a horizontal line to divide the plot by gene
                             "9" = gpar(lty = 2)), # a horizontal line to divide the plot by gene
           txt_gp = fpTxtGp(summary = list(gpar(cex=1, fontface = "plain")), # text size etc
                            label = gpar(cex=1),
                            ticks = gpar(cex=1),
                            xlab  = gpar(fontface="bold", cex = 1)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0.15, ...)} # format the "summary" rows to suit the main rows
)
dev.off()
