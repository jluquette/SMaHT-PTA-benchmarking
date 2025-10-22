#!/bin/bash

batch_list=("batch-1_lung" "batch-2_colon")

lambda_list=("1024 512" "128 256" "64 32")

for batch in "${batch_list[@]}"; do
    for lambda in "${lambda_list[@]}"; do
        for multisample in "true" "false"; do
            for binsize in "500000" "100000"; do
                echo "${batch} submitted with ${lambda} with ${multisample} binsize ${binsize}"
                sbatch sbatch_submit_binsize_one-lambda_hiscanner.sh $batch $multisample $binsize "ST002" "${lambda}"
            done
        done
    done
done

batch="batch-3_scan2-neuron"
lambda_list=("1024 512 128 256 64 32")
echo "Calling HiScanner on SCAN2 neuron (control)"
while read sample; do
    for lambda in "${lambda_list[@]}"; do
        for multisample in "true" "false"; do
            for binsize in "500000" "100000"; do
                echo "${batch} submitted with ${lambda} with ${multisample} binsize ${binsize} ($sample)"
                sbatch sbatch_submit_binsize_one-lambda_hiscanner.sh $batch $multisample $binsize $sample "${lambda}"
            done
        done
    done
done < "../batch-3_scan2-neuron/samples.txt"
