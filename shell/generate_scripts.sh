
########## GENERATE SCRIPTS THAT CAN BE RUN IN PARALLEL ##########

# generate generic script for MR – clumping at r2 < 0.2
echo "
EXPOSURE_DATA<-'${EXPOSURE_DATA}'
OUTCOME<-'${OUTCOME}'
" > scripts/script_${EXPOSURE_DATA}_${OUTCOME}.R
            
cat read_exposure_data_${EXPOSURE_DATA}.R >> scripts/script_${EXPOSURE_DATA}_${OUTCOME}.R

cat read_outcome_data_${OUTCOME}.R >> scripts/script_${EXPOSURE_DATA}_${OUTCOME}.R

cat mr_analysis.R >> scripts/script_${EXPOSURE_DATA}_${OUTCOME}.R
