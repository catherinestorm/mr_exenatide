setwd('~/Documents/UCL_PhD/PhD_Project/mr_exenatide_pd/version_2020')

EXPOSURE_DATA<-'eqtlgen'
OUTCOME<-'t2dm_risk'


# read in packages

library("TwoSampleMR")
library("dplyr")
library("stringr")
library("readr")
source("~/Documents/UCL_PhD/PhD_Project/mr_resources/papers_resources/methods/mr_fclr/FCLR.R")

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

dat <- clump_data(dat, clump_r2 = 0.4)

dat2 <- dat_to_MRInput(dat, get_correlation=TRUE)


ivw <- MendelianRandomization::mr_ivw(dat2[[1]], correl=TRUE)

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
                 
                 r = 5)

ivw
fclr.res









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

pca <- data.frame("exposure_data" = EXPOSURE_DATA, "outcome" = OUTCOME, "ensembl" = "ENSG00000112164", "beta" = beta_IVWcorrel.pc, "se" = se_IVWcorrel.fixed.pc, "p" = pvalue, "het" = pvalue.heter.stat)

pca









# compare to classical methods
# ENSG00000112164


# liberal analysis

dat0 <- harmonise_data(
  exposure_dat = exp,
  outcome_dat = out)

dat <- subset(dat0, dat0$mr_keep == TRUE)


dat0.2 <- clump_data(dat, clump_r2 = 0.6)


#STEIGER TEST

steiger_nalls0.2 <- directionality_test(dat0.2)

steiger_nalls0.2_1 <- subset(steiger_nalls0.2, correct_causal_direction == TRUE)

dat_steiger0.2 <- subset(dat0.2, dat0.2$exposure %in% steiger_nalls0.2_1$exposure)





# correlated where > 2 SNPs available

dat_steiger_keep <- dat_steiger0.2


  
  unique_exposures0 <- unique(dat_steiger_keep$exposure)
  
  unique_exposures0 <- sort(unique_exposures0)
  
  unique_exposures <- unique_exposures0[!is.na(unique_exposures0)]
  
  
  corr_results <- data.frame()
  
  corr_egger_intercept <- data.frame()
  
  
  
  for (i in 1:length(unique_exposures)) {
    dat_steiger_keep1 <- dat_steiger_keep[which(dat_steiger_keep$exposure == unique_exposures[i]), ] #subset data to keep only 1 exposure of interest
    
    dat2 <- dat_to_MRInput(dat_steiger_keep1, get_correlation=TRUE)
    
    ivw <- MendelianRandomization::mr_ivw(dat2[[1]], correl=TRUE)
    egger <- MendelianRandomization::mr_egger(dat2[[1]], correl=TRUE)
    maxlik <- MendelianRandomization::mr_maxlik(dat2[[1]], correl=TRUE)
    temp_data <- rbind(c(ivw@Outcome, ivw@SNPs, ivw@class[1], ivw@Estimate, ivw@StdError, ivw@Pvalue),
                       c(egger@Outcome, egger@SNPs, egger@class[1], egger@Estimate, egger@StdError.Est, egger@Pvalue.Est),
                       c(maxlik@Outcome, maxlik@SNPs, maxlik@class[1], maxlik@Estimate, maxlik@StdError, maxlik@Pvalue))
    
    corr_results <- rbind(corr_results, temp_data)
    
    temp_data1 <- cbind(egger@Outcome, egger@Intercept, egger@CILower.Int, egger@CIUpper.Int, egger@Pvalue.Int, egger@Heter.Stat[1], egger@Heter.Stat[2])
    corr_egger_intercept <- rbind(corr_egger_intercept, temp_data1)
    
    
  }
  
  
  
  #correlated results
  corr_results1 <- corr_results
  names(corr_results1) <- c("outcome", "nsnp", "method", "beta", "se", "p")
  
  corr_results1[which(corr_results1$method == "IVW"), "exposure"] <- unique_exposures
  corr_results1[which(corr_results1$method == "Egger"), "exposure"] <- unique_exposures
  corr_results1[which(corr_results1$method == "MaxLik"), "exposure"] <- unique_exposures
  corr_results1[,4:6] <- sapply(corr_results1[,4:6], function(x){as.numeric(as.character(x))})
  
  names(corr_egger_intercept) <- c("outcome","egger_intercept","lower_ci", "upper_ci", "pvalue", "Q", "Q_pval")
  
  corr_results1
  corr_egger_intercept
  
  
  
  
  
  
  
  
  
  
  # T2DM
  
  # 0.2 corr_results1
  outcome nsnp method        beta         se          p        exposure
  1 T2DM risk    7    IVW -0.19477704 0.07850813 0.01310220 ENSG00000112164
  2 T2DM risk    7  Egger -0.03858583 0.16234524 0.81213114 ENSG00000112164
  3 T2DM risk    7 MaxLik -0.20303027 0.08170154 0.01295428 ENSG00000112164
  
  
  # 0.4 corr_results1
  outcome nsnp method       beta         se           p        exposure
  1 T2DM risk    9    IVW -0.2118562 0.07831852 0.006829243 ENSG00000112164
  2 T2DM risk    9  Egger -0.0436827 0.14601045 0.764806425 ENSG00000112164
  3 T2DM risk    9 MaxLik -0.2301759 0.08237103 0.005199937 ENSG00000112164
  
  
  # 0.6 corr_results1
  outcome nsnp method        beta         se           p        exposure
  1 T2DM risk   14    IVW -0.20191211 0.08513517 0.017708135 ENSG00000112164
  2 T2DM risk   14  Egger -0.08670134 0.10292077 0.399559307 ENSG00000112164
  3 T2DM risk   14 MaxLik -0.24716729 0.09263313 0.007625007 ENSG00000112164
  
  
  
  
  
  # FCLR r = 2; no clump
  # ~ 202 SNPs
  $fliml.est
  [1] -0.3217164
  
  $fliml.se
  [1] 0.1229496
  
  $AR.p
  [1] 0.01232414
  
  $LM.p
  [1] 0.006367888
  
  $CLR.p
  [1] 0.006255805
  
  $CLR.CI
  [1] -0.60 -0.09
  
  
  # FCLR r = 4; clump @ 0.6
  # ~ 14 SNPs
  $fliml.est
  [1] -0.2500722
  
  $fliml.se
  [1] 0.08880558
  
  $AR.p
  [1] 0.025526
  
  $LM.p
  [1] 0.004422893
  
  $CLR.p
  [1] 0.004442475
  
  $CLR.CI
  [1] -0.44 -0.08
  

  
  
  # FCLR r = 4; clump @ 0.8 
  # ~ 20 SNPs
  $fliml.est
  [1] -0.2119091
  
  $fliml.se
  [1] 0.08536469
  
  $AR.p
  [1] 0.02326305
  
  $LM.p
  [1] 0.01457817
  
  $CLR.p
  [1] 0.01373006
  
  $CLR.CI
  [1] -0.39 -0.05
  
  
  
  
  
  
  # PCA no clump
  exposure_data   outcome         ensembl       beta        se           p        het
  1       eqtlgen t2dm_risk ENSG00000112164 -0.2774226 0.1015626 0.006303836 0.0726697
  
  
  # PCA clump @ 0.6
  exposure_data   outcome         ensembl      beta         se           p       het
  1       eqtlgen t2dm_risk ENSG00000112164 -0.207165 0.07923327 0.008932647 0.2131153
  
  
  # PCA clump @ 0.8
  exposure_data   outcome         ensembl       beta         se            p          het
  1       eqtlgen t2dm_risk ENSG00000112164 -0.3206529 0.06608039 1.219357e-06 0.0004947439
  
  
  
  
  
  
