

EXPOSURE_DATA<-'eqtlgen'
OUTCOME<-'t2dm_risk'


# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")
library("readr")
source("FCLR_with_pvalue.R")
library("forestplot")


# download LD reference panel for each of the 5 super-populations in the 1000 genomes reference dataset. e.g. for the European super population it has the following files:
# from https://mrcieu.github.io/ieugwasr/articles/local_ld.html
devtools::install_github("explodecomputer/genetics.binaRies")
genetics.binaRies::get_plink_binary()

# calculated maf frq in plink 
# ./plink --bfile EUR --freq --out EUR
# cat EUR.frq | sed -r 's/^\s+//g' | sed -r 's/\s+/,/g' > EUR.frq.csv

eur_maf <- read_csv("EUR.frq.csv")


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

# keep GLP1R only
# GLP1R == ENSG00000112164
exp <- subset(exp0, (exp0$exposure == "ENSG00000112164"))



### load outcome data

out <- read_outcome_data(snps = exp$SNP,
                         filename = "data/xue2018_t2dm_risk_with_eqtl_rsid.txt",
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

# clump
clump_thresh <- 0.6
dat <- clump_data(dat, clump_r2 = clump_thresh)

# get LD matrix
dat2 <- dat_to_MRInput(dat, get_correlation=TRUE)

# IVW correcting for LD
ivw <- MendelianRandomization::mr_ivw(dat2[[1]], correl=TRUE)

# plot
mrbase<- mr(dat, method_list = "mr_ivw")
p1 <- mr_scatter_plot(mrbase, dat)
for (i in 1:length(p1)){
  ggplot2::ggsave(p1[[i]], file="results/plots/scat_glp1r_t2dm_fclr_pca_ivw_r2_0.6.pdf", width=7, height=7)
}




# FCLR and FLIML method

# get the maf in europeans for the SNPs we have
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





# scree plot to decide r
pcs <- prcomp(dat2[[1]]@correlation)
screeplot(pcs,
          type = "lines")



# fclr and LIML
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


# put ivw, fclr and fliml into one table
ivw_for_table <- data.frame("clump_thresh" = clump_thresh, "nsnp" = ivw@SNPs, "method"=ivw@class[1], "beta"=ivw@Estimate, "se"=ivw@StdError, "p"=ivw@Pvalue, "lo_ci" = (ivw@Estimate-1.96*ivw@StdError),  "up_ci" = (ivw@Estimate+1.96*ivw@StdError),  "het" = NA)
fclr.res_for_table1 <- data.frame("clump_thresh" = clump_thresh,"nsnp" = NA, "method"="F-LIML", "beta"=fclr.res$fliml.est[1], "se"=fclr.res$fliml.se[1], "p"=fclr.res$fliml.pval[1], "lo_ci" = (fclr.res$fliml.est[1]-1.96*fclr.res$fliml.se[1]),  "up_ci" = (fclr.res$fliml.est[1]+1.96*fclr.res$fliml.se[1]),   "het" = NA)
fclr.res_for_table2 <- data.frame("clump_thresh" = clump_thresh,"nsnp" = NA, "method"="F-CLR", "beta"=NA, "se"=NA, "p"=fclr.res$CLR.p[1], "lo_ci" = fclr.res$CLR.CI[1],  "up_ci" = fclr.res$CLR.CI[2],   "het" = NA)
fclr.res_for_table <- rbind(fclr.res_for_table1,fclr.res_for_table2)

temp <- rbind(ivw_for_table,fclr.res_for_table)






# IVWPC method
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

pca <- data.frame("clump_thresh" = clump_thresh,"nsnp" = NA, "method"="IVWPC", "beta"=beta_IVWcorrel.pc, "se"=se_IVWcorrel.fixed.pc, "p"=pvalue, "lo_ci" = (beta_IVWcorrel.pc-1.96*se_IVWcorrel.fixed.pc),  "up_ci" = (beta_IVWcorrel.pc+1.96*se_IVWcorrel.fixed.pc), "het" = pvalue.heter.stat)

full_data <- rbind(temp,pca)

full_data$or <- exp(full_data$beta)
full_data$or_lci95 <- exp(full_data$lo_ci)
full_data$or_uci95 <- exp(full_data$up_ci)


write.table(full_data, "results/full_res_r2_0.6_fclr_pca_ivw.txt", row.names = F, sep = "\t")





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

# format the data for the forest plot
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

# DATA TO DISPLAY AS NUMBER IN THE FOREST PLOT
table_text <- data[, c("method","nsnp", "or_ci", "roundp")]
table_text <- rbind(c("Method","No. SNPs", "OR (95% CI)", "FDR-corrected P"), table_text)




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


