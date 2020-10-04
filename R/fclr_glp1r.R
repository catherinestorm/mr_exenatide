setwd('~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020')

EXPOSURE_DATA<-'eqtlgen'
OUTCOME<-'t2dm_risk'


# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")
library("readr")
#source("~/Documents/UCL_PhD/PhD_Project/mr_resources/papers_resources/methods/mr_fclr/FCLR.R")
source("~/Documents/UCL_PhD/PhD_Project/mr_resources/papers_resources/methods/mr_fclr/FCLR_with_pvalue.R")


eur_maf <- read_csv("~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020/for_clumping/EUR.frq.csv")


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

exp0 <- subset(exp0, exp0$pval.exposure < 5e-5)

#exp_to_keep <- read_csv("genes_for_analysis.csv")
exp <- subset(exp0, (exp0$exposure == "ENSG00000112164"))



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


# harmonise

dat0 <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out)




dat <- subset(dat0, dat0$mr_keep == TRUE)

clump_thresh <- 0.6
dat <- clump_data(dat, clump_r2 = clump_thresh)

dat2 <- dat_to_MRInput(dat, get_correlation=TRUE)


ivw <- MendelianRandomization::mr_ivw(dat2[[1]], correl=TRUE)

mrbase<- mr(dat, method_list = "mr_ivw")
p1 <- mr_scatter_plot(mrbase, dat)
for (i in 1:length(p1)){
  ggplot2::ggsave(p1[[i]], file="results/plots/scat_glp1r_t2dm_fclr_pca_ivw_r2_0.6.pdf", width=7, height=7)
}


#ld_matr<-ld_matrix(dat$SNP)

#dat2 <- dat_to_MRInput(dat[!(dat$SNP == "rs1678690"),], get_correlation=TRUE)


# calculated frq in plink on google cloud engine pdtreatment
# ./plink --bfile for_clumping/EUR --freq --out for_clumping/EUR
# cat for_clumping/EUR.frq | sed -r 's/^\s+//g' | sed -r 's/\s+/,/g' > for_clumping/EUR.frq.csv
# cd /Users/catherinestorm/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020/for_clumping
# gcloud compute scp catherinestorm@pdtreatment://home/catherinestorm/for_clumping/EUR.frq.csv .




eur_maf_keep0 <- subset(eur_maf, eur_maf$SNP %in% dat2[[1]]@snps)
eur_maf_keep <- as.data.frame(eur_maf_keep0)

length(dat2[[1]]@snps) == nrow(eur_maf_keep)
all(dat2[[1]]@snps == eur_maf_keep$SNP)

rownames(eur_maf_keep) <- eur_maf_keep$SNP

eur_maf_keep <- eur_maf_keep[dat2[[1]]@snps,]

all(dat2[[1]]@snps == eur_maf_keep$SNP)



mafs <- dat[,c("SNP","eaf.exposure","eaf.outcome")]
row.names(mafs) <- mafs$SNP
mafs <- mafs[dat2[[1]]@snps,]
all(dat2[[1]]@snps == mafs$SNP)





# scree plot
pcs <- prcomp(dat2[[1]]@correlation)
screeplot(pcs,
          type = "lines")




fclr.res <- FCLR(nx = dat$samplesize.exposure,
                 ny = dat$samplesize.outcome,
                 
                 bx = dat2[[1]]@betaX,
                 se.bx = dat2[[1]]@betaXse,
                 
                 by = dat2[[1]]@betaY,
                 se.by = dat2[[1]]@betaYse,
                 
                 mafs.x = mafs$eaf.exposure,
                 mafs.y = mafs$eaf.outcome,
                 mafs.z = eur_maf_keep$MAF,
                 
                 LD=dat2[[1]]@correlation,
                 
                 r = 3)

ivw_for_table <- data.frame("clump_thresh" = clump_thresh, "nsnp" = ivw@SNPs, "method"=ivw@class[1], "beta"=ivw@Estimate, "se"=ivw@StdError, "p"=ivw@Pvalue, "lo_ci" = (ivw@Estimate-1.96*ivw@StdError),  "up_ci" = (ivw@Estimate+1.96*ivw@StdError),  "het" = NA)
fclr.res_for_table1 <- data.frame("clump_thresh" = clump_thresh,"nsnp" = NA, "method"="F-LIML", "beta"=fclr.res$fliml.est[1], "se"=fclr.res$fliml.se[1], "p"=fclr.res$fliml.pval[1], "lo_ci" = (fclr.res$fliml.est[1]-1.96*fclr.res$fliml.se[1]),  "up_ci" = (fclr.res$fliml.est[1]+1.96*fclr.res$fliml.se[1]),   "het" = NA)
fclr.res_for_table2 <- data.frame("clump_thresh" = clump_thresh,"nsnp" = NA, "method"="F-CLR", "beta"=NA, "se"=NA, "p"=fclr.res$CLR.p[1], "lo_ci" = fclr.res$CLR.CI[1],  "up_ci" = fclr.res$CLR.CI[2],   "het" = NA)
fclr.res_for_table <- rbind(fclr.res_for_table1,fclr.res_for_table2)

temp <- rbind(ivw_for_table,fclr.res_for_table)






# PCA
rho<-ld_matrix(dat$SNP)

dat_keep <- subset(dat, dat$SNP %in% sub("\\_.*","", colnames(rho))) #since some SNPs may have been lost in ld matrix generation


betaXG =dat_keep$beta.exposure
sebetaXG=dat_keep$se.exposure
betaYG=dat_keep$beta.outcome
sebetaYG =dat_keep$se.outcome

Phi = (betaXG/sebetaYG)%o%(betaXG/sebetaYG)*rho 
#summary(prcomp(Phi, scale=FALSE)) 
K = which(cumsum(prcomp(Phi, scale=FALSE)$sdev^2/sum((prcomp(Phi, scale=FALSE)$sdev^2)))>0.99)[1]
# K is number of principal components to include in analysis 
# this code includes principal components to explain 99% of variance in the risk factor 

betaXG0 = as.numeric(betaXG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]) 
betaYG0 = as.numeric(betaYG%*%prcomp(Phi, scale=FALSE)$rotation[,1:K]) 
Omega = sebetaYG%o%sebetaYG*rho 
pcOmega = t(prcomp(Phi, scale=FALSE)$rotation[,1:K])%*%Omega%*%prcomp(Phi, scale=FALSE)$rotation[,1:K] 
beta_IVWcorrel.pc = solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0)*t(betaXG0)%*%solve(pcOmega)%*%betaYG0
se_IVWcorrel.fixed.pc = sqrt(solve(t(betaXG0)%*%solve(pcOmega)%*%betaXG0))
pvalue <- 2 * pnorm(-abs(beta_IVWcorrel.pc/se_IVWcorrel.fixed.pc))
#pvalue <- 2 * pt(-abs(beta_IVWcorrel.pc/se_IVWcorrel.fixed.pc), df = length(betaXG) -1)#t-dist

#Q statistic (accounting for correlation) using principal components: 
rse = betaYG0- beta_IVWcorrel.pc*betaXG0 
rse.corr = sqrt(t(rse)%*%solve(pcOmega)%*%rse/(K-1)) # 
heter.stat <- (K - 1)*(rse.corr^2)
pvalue.heter.stat <- pchisq(heter.stat, df = K-1, lower.tail = FALSE)

pca <- data.frame("clump_thresh" = clump_thresh,"nsnp" = NA, "method"="PCA", "beta"=beta_IVWcorrel.pc, "se"=se_IVWcorrel.fixed.pc, "p"=pvalue, "lo_ci" = (beta_IVWcorrel.pc-1.96*se_IVWcorrel.fixed.pc),  "up_ci" = (beta_IVWcorrel.pc+1.96*se_IVWcorrel.fixed.pc), "het" = pvalue.heter.stat)

full_data <- rbind(temp,pca)


write.table(full_data, "results/full_res_r2_0.6_fclr_pca_ivw.txt", row.names = F, sep = "\t")


full_data$or <- exp(full_data$beta)
full_data$or_lci95 <- exp(full_data$lo_ci)
full_data$or_uci95 <- exp(full_data$up_ci)


# FOREST PLOT
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


data <- full_data

data$or_ci <- str_c(format_numbers(data$or), " (", format_numbers(data$or_lci95), ", ",format_numbers(data$or_uci95), ")")
data$roundp <- lapply(X = data$p, FUN = format_numbers)

data$method <- as.character(data$method)

# dummy OR
data$or[data$method == "F-CLR"] <- data$or_lci95[data$method == "F-CLR"]
data$or_ci[data$method == "F-CLR"] <- str_c("NA (", format_numbers(data$or_lci95), ", ",format_numbers(data$or_uci95), ")")


# DATA FOR THE CREATING THE FOREST PLOT
forest_data <- data[,c("or", "or_lci95", "or_uci95", "method")]
forest_data <- forest_data[complete.cases(forest_data),]
forest_data <- rbind(c(NA, NA, NA, NA),forest_data) # size we have a header row we want a blank ros at the top

# DATA TO DISPLAY AS NUMBETS IM THE THE FOREST PLOT
table_text <- data[, c("method","nsnp", "or_ci", "roundp")]
table_text <- rbind(c("Method","No. SNPs", "OR (95% CI)", "FDR-adjusted P"), table_text)




pdf("results/plots/forest_glp1r_t2dm_fclr_pca_ivw_r2_0.6.pdf", width = 23,height = 10, onefile=FALSE)
forestplot(table_text, # columns to include
           graph.pos = 3, # where graph is
           boxsize = 0.1, # size of the box
           align = "l", # align text to left 
           mean = forest_data$or,
           lower = forest_data$or_lci95,
           upper = forest_data$or_uci95,
           xlog=TRUE, # log scale if odds ratios
           is.summary= c(forest_data$method %in% c("F-CLR", "NA")),
           col=fpColors(box="darkred", line = "darkred", summary = "darkred"), # colours for blood vs brain tissue
           xlab="T2DM odds ratio",
           hrzl_lines = list("2" = gpar(lty = 2)), # a horizontal line at the top under the header
           txt_gp = fpTxtGp(summary = list(gpar(cex=3, fontface = "plain")),
                            label = gpar(cex=3),
                            ticks = gpar(cex=3),
                            xlab  = gpar(fontface="bold", cex = 3)),
           fn.ci_sum=function(col, size, ...) {fpDrawNormalCI(clr.line = col, clr.marker = col, size=0, ...)} # format the "summary" rows to suit the main rows
)
dev.off()





# 5e-5
clump_thresh nsnp method       beta         se           p       clr_ci       clr_p       het
1          0.6   14    IVW -0.2019121 0.08513517 0.017708135         <NA>          NA        NA
2          0.6   NA F-LIML -0.2478864 0.08903123          NA -0.43, -0.08 0.004748768        NA
3          0.6   NA    PCA -0.2071650 0.07923327 0.008932647         <NA>          NA 0.2131153



# 5e-8
clump_thresh nsnp method       beta         se          p      clr_ci       clr_p        het
1          0.6    4    IVW -0.1879417 0.10955510 0.08625370        <NA>          NA         NA
2          0.6   NA F-LIML -0.3394059 0.11916541         NA -0.6, -0.12 0.002765291         NA
3          0.6   NA    PCA -0.1914234 0.08060622 0.01755862        <NA>          NA 0.07357023

  
  
  
