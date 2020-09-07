echo -n > failed_scripts_liberal_logs.txt

while read EXPOSURE_DATA; do
    while read OUTCOME; do

        grep -L "mission_complete" scripts/nohup_script_${EXPOSURE_DATA}_${OUTCOME}.log >> failed_scripts_liberal_logs.txt


    done < outcomes.txt
done < exposure_data.txt

cat failed_scripts_liberal_logs.txt | sed "s/nohup_//g" | sed "s/.log//g" > failed_scripts_liberal.txt


while read FAILED; do
    export LOG=$(sed "s/script_/nohup_script_/g" <<< ${FAILED})
    nohup Rscript ${FAILED}.R &> ${LOG}.log &
    wait
done < failed_scripts_liberal.txt
