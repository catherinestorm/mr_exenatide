# mr_exenatide



## Introduction
These scripts were used in a Mendelian randomization analysis of exenatide in Parkinson's disease. Full methods can be found [here]().

## Citation
If you use the code, please cite: [Storm CS, Kia DA, Almramhi M, Athauda D, International Parkinson’s Disease Genomics Consortium (IPDGC), Burgess S, Foltynie T, Wood NW. "Predicting the efficacy of exenatide in Parkinson’s disease using genetics – a Mendelian randomization study." bioRxiv. 2020.]()




## Pipeline Overview



1. Set up work space. Specify exposure and outcome datasets.

```mkdir data
mkdir data
mkdir data/progression
mkdir results
mkdir results/plots
mkdir scripts

echo "eqtlgen
psychencode" > exposure_data.txt


echo "t2dm_risk
bmi" > non_pd_outcomes.txt

echo "pd_risk
pd_aao" > pd_outcomes.txt

echo "cont_HY
cont_MMSE
cont_MOCA
cont_SEADL
cont_UPDRS1_scaled
cont_UPDRS2_scaled
cont_UPDRS3_scaled
cont_UPDRS4_scaled
cont_UPDRS_scaled
surv_DEMENTIA
surv_DEPR
surv_DYSKINESIAS
surv_HY3" > progression_outcomes.txt

cat progression_outcomes.txt >> pd_outcomes.txt


while read EXPOSURE_DATA; do
    while read OUTCOME; do
        mkdir results/plots/${EXPOSURE_DATA}_${OUTCOME}
    done < pd_outcomes.txt
done < exposure_data.txt


while read OUTCOME; do
    mkdir results/plots/eqtlgen_${OUTCOME}
done < non_pd_outcomes.txt


```


2. Download data. Prepare the data for the Mendelian randomization analysis.
To download the Parkinson's disease risk, age at onset and progression GWAS data, go to the below links:

Download Parkinson's disease risk GWAS data: https://bit.ly/2ofzGrk
Download Parkinson's disease age at onset GWAS data: http://pdgenetics.org/resources
Download Parkinson's disease progression GWAS data: https://pdgenetics.shinyapps.io/pdprogmetagwasbrowser/
NOTE: make sure to untick the P < 0.05 box in order to get the full summary statistics for progression

```
nohup bash ./mr_exenatide/shell/data_download.sh &> nohup_data_download.log &

nohup Rscript ./mr_exenatide/R/data_prep_eqtl.R &> nohup_data_prep_eqtl.log &

nohup Rscript ./mr_exenatide/R/data_prep_outcome_gwas.R &> nohup_data_prep_outcome_gwas.log &

```

3. Generate scripts to run the analysis.

```
# generate scripts for non-Parkinson's outcomes
while read OUTCOME; do
    export EXPOSURE_DATA="eqtlgen"
    export OUTCOME=${OUTCOME}
    bash ./mr_exenatide/shell/generate_scripts.sh
done < non_pd_outcomes.txt


# generate read_outcome_data scripts for progression
while read OUTCOME; do
    cat ./mr_exenatide/R/read_outcome_data_progression.R > read_outcome_data_${OUTCOME}.R
done < progression_outcomes.txt


# generate scripts for Parkinson's-related outcomes
while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        bash ./mr_exenatide/shell/generate_scripts.sh
    done < pd_outcomes.txt
done < exposure_data.txt

```

4. Run the scripts
```
nohup bash ./mr_exenatide/shell/run_scripts.sh &> nohup_run_scripts.log &

nohup Rscript ./mr_exenatide/R/mr_ivw_fclr_fliml_ivwpc_glp1r_t2dm.R &> nohup_mr_ivw_fclr_fliml_ivwpc_glp1r_t2dm.log &

```

5. Combine results from all analyses
```
Rscript ./mr_exenatide/R/combine_data.R
```

6. Generate forest plots.
```
Rscript ./mr_exenatide/R/make_forest_plots.R
```
