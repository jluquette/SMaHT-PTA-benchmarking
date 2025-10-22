#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH -p priopark
#SBATCH -A park_contrib


# Step 1 - generate the unfiltered, primary tables of mutations
if [ ! -f all_SCAN2_pass_snv_indel_dnv.csv.gz ]; then
    echo "Generating all_SCAN2_pass_snv_indel_dnv.csv.gz..."
    Rscript -e 'suppressMessages(library(scan2)); setwd("objs/summary"); s <- load.summary(); fwrite(passing(s, combine.mnv=TRUE)[muttype != "mnv"], file="../../all_SCAN2_pass_snv_indel_dnv.csv")'
    gzip all_SCAN2_pass_snv_indel_dnv.csv
else
    echo "Skipping all_SCAN2_pass_snv_indel_dnv.csv.gz..."
fi

if [ ! -f all_SCAN2_pass_snv_indel_dnv.relaxed.csv.gz ]; then
    echo "Generating all_SCAN2_pass_snv_indel_dnv.relaxed.csv.gz..."
    Rscript -e 'suppressMessages(library(scan2)); setwd("objs/summary_bulk_binom_prob_1e-5"); s <- load.summary(); fwrite(passing(s, combine.mnv=TRUE)[muttype != "mnv"], file="../../all_SCAN2_pass_snv_indel_dnv.relaxed.csv")'
    gzip all_SCAN2_pass_snv_indel_dnv.relaxed.csv
else
    echo "Skipping all_SCAN2_pass_snv_indel_dnv.relaxed.csv.gz..."
fi
