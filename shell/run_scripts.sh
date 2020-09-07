while read EXPOSURE_DATA; do
    while read OUTCOME; do
        export EXPOSURE_DATA=${EXPOSURE_DATA}
        export OUTCOME=${OUTCOME}
        nohup Rscript scripts/script_${EXPOSURE_DATA}_${OUTCOME}.R &> scripts/nohup_script_${EXPOSURE_DATA}_${OUTCOME}.log &
        wait
    done < outcomes.txt
done < exposure_data.txt
