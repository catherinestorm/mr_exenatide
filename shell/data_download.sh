# download eQTLGen Data

wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz -P data

wget https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz -P data



# download PsychENCODE Data
wget http://resource.psychencode.org/Datasets/Derived/QTLs/SNP_Information_Table_with_Alleles.txt -P data

wget http://resource.psychencode.org/Datasets/Derived/QTLs/DER-08a_hg19_eQTL.significant.txt -P data



# download body mass index GWAS data
wget https://portals.broadinstitute.org/collaboration/giant/images/c/c8/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz -P data

gunzip Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz



# download type 2 diabetes mellitus GWAS data
wget https://cnsgenomics.com/data/t2d/Xue_et_al_T2D_META_Nat_Commun_2018.gz
gunzip Xue_et_al_T2D_META_Nat_Commun_2018.gz



# download Parkinson's disease risk GWAS data
# go to https://bit.ly/2ofzGrk and download the file
unzip nallsEtAl2019_excluding23andMe_allVariants.tab.zip


# download Parkinson's disease age at onset GWAS data
# go to https://bit.ly/2ofzGrk and download the relevant data
# http://pdgenetics.org/resources



# download Parkinson's disease progression GWAS data
# go to the below link and download the relevant data
# MAKE SURE TO UNTICK WITH P < 0.05 BOX IN ORDER TO GET THE FULL SUMMARY STATISTICS
# https://pdgenetics.shinyapps.io/pdprogmetagwasbrowser/

