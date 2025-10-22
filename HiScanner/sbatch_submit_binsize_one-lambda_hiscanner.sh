#!/bin/bash
#SBATCH -n 1                 
#SBATCH -t 5-00:00:00
#SBATCH -p park,medium
#SBATCH --mem=16G             
#SBATCH --account=park
#SBATCH -o logs/cnv_hiscanner_one-lambda_%j.out
#SBATCH -e logs/cnv_hiscanner_one-lambda_%j.err

# Load modules
module load gcc/14.2.0
module load conda/miniforge3/24.11.3-0
module load samtools

batch_name=$1
multisample=$2
binsize=$3 #500000
sample=$4 #"ST002"

if [ "$multisample" = "true" ]; then
    multisample_suffix="_multisample-seg"
else
    multisample_suffix=""
fi

source setup_variables.sh ${batch_name}

WORKING_DIR=${PTA_DIR}/${batch_name}


BINSIZE=$binsize
LAMBDA_LIST=($5) #"1 2 4" #(2 4 8 16 32 64 128 256 512 1024)

lambda_str=$(IFS=- ; echo "${LAMBDA_LIST[*]}")

HISCANNER_DIR_BASE=${WORKING_DIR}/hiscanner_multisample-${multisample}
HISCANNER_DIR_BASE=${HISCANNER_DIR_BASE}/binsize_${binsize}_lambda_${lambda_str}
mkdir -p $HISCANNER_DIR_BASE
HISCANNER_DIR=${HISCANNER_DIR_BASE}/${sample}
mkdir -p $HISCANNER_DIR

echo $HISCANNER_DIR

conda run -n hiscanner_test hiscanner init --output $HISCANNER_DIR


output_file=${HISCANNER_DIR}/"metadata.txt"
echo -e "bamID\tbam\tsinglecell" > "$output_file"

 # Process bulk samples
for bulk_bam in "${WORKING_DIR}/bulk_bams/${sample}/"*.bam; do
    if [ -f "$bulk_bam" ]; then
        # Extract sample name from BAM header using samtools
        sample_name=$(samtools view -H "$bulk_bam" | grep '@RG' | grep -o 'SM:[^[:space:]]*' | head -1 | cut -d: -f2)

        # If samtools couldn't extract a name, use the filename as fallback
        if [ -z "$sample_name" ]; then
            sample_name=$(basename "$bulk_bam" .bam)
        fi

        # Add to metadata file with singlecell=N
        echo -e "${sample_name}\t${bulk_bam}\tN" >> "$output_file"
    fi
done


for sc_bam in "${WORKING_DIR}/sc_bams/${sample}/"*.bam; do
    if [ -f "$sc_bam" ]; then
        # Extract sample name from BAM header using samtools
        cell_id=$(samtools view -H "$sc_bam" | grep '^@RG' | sed -n 's/.*\tSM:\([^ \t]*\).*/\1/p' | head -1)
        # If samtools couldn't extract a name, create one based on filename and sample
        if [ -z "$cell_id" ]; then
            cell_id="${sample}_$(basename "$sc_bam" .bam)"
        fi

        echo -e "${cell_id}\t${sc_bam}\tY" >> "$output_file"
    fi
done

            
cat <<EOT > "${HISCANNER_DIR}/config_orig.yaml"
scan2_output: ${WORKING_DIR}/scan2/${sample}
metadata_path: ${HISCANNER_DIR}/metadata.txt
outdir: ${HISCANNER_DIR}/hiscanner_output
use_multisample_segmentation: ${multisample}

# Reference data and external tools - MUST BE SET
fasta_folder: /n/data1/hms/dbmi/park/ann_caplin/HiScanner/GRCh38/split
mappability_folder_stem: /n/data1/hms/dbmi/park/ann_caplin/HiScanner/GRCh38/hg38_mappability/150mer.

# Analysis parameters - defaults provided but can be modified
rdr_only: false
binsize: ${BINSIZE}
chrom_list: ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
lambda_range: [$(IFS=,; echo "${LAMBDA_LIST[*]}" | sed 's/,/, /g')]
lambda_value: ${LAMBDA_LIST}
max_wgd: 1
batch_size: 5            # Number of cells to process in parallel
depth_filter: 0          # Minimum read depth filter
ado_threshold: 0.2       # Allelic dropout threshold
threads: 4               # Number of parallel threads to use
ado_plot_baf_distribution: true  # Set to false to disable BAF distribution plots

# ADO analysis options
aggregate_every_k_snp: true  # Whether to aggregate every k SNPs for ADO analysis
k: 5

# Debug and performance options
keep_raw_files: false    # Whether to keep intermediate raw files (default: false)
rerun: false            # Force rerun of existing results (default: false)
keep_temp_files: false  # Keep temporary files from multisample segmentation


# Cluster options
use_cluster: false      # Whether to use cluster for computation
cluster_queue: "park"   # Cluster queue to use
max_cluster_jobs: 100   # Maximum number of concurrent cluster jobs

EOT

cp ${HISCANNER_DIR}/config_orig.yaml ${HISCANNER_DIR}/config.yaml

cat <<EOT > "${HISCANNER_DIR}/cluster_config_orig.yaml"
__default__:
  account: "park"
  partition: "park,short"
  time: "12:00:00"  # HH:MM:SS format
  nodes: 1
  ntasks: 1
  mem: "4G"
  job-name: "hiscanner_{rule}"
  output: "cluster_logs/{rule}-%j.out"
  error: "cluster_logs/{rule}-%j.err"

get_read_pos:
  # Typically lighter on resources but IO-intensive
  time: "12:00:00"
  mem: "4G"
  ntasks: 1
  partition: "park,short"

run_bicseq_norm:
  # More computationally intensive
  time: "12:00:00"
  mem: "32G"
  ntasks: 8
  partition: "park,short"
  
EOT


cp ${HISCANNER_DIR}/cluster_config_orig.yaml ${HISCANNER_DIR}/cluster.yaml

cd $HISCANNER_DIR

conda run -n hiscanner_test hiscanner run --step snp      # Check SCAN2 results
conda run -n hiscanner_test hiscanner run --step phase    # Process SCAN2 results
conda run -n hiscanner_test hiscanner run --step ado      # ADO analysis to identify optimal bin size

conda run -n hiscanner_test hiscanner --config config.yaml run --step normalize --use-cluster
conda run -n hiscanner_test hiscanner run --step segment  # Segmentation

conda run -n hiscanner_test hiscanner run --step cnv



mv ${HISCANNER_DIR}/hiscanner_output/final_calls ${HISCANNER_DIR}/hiscanner_output/final_calls_old

for lambda in "${LAMBDA_LIST[@]}"; do
    if ! [ -d ${HISCANNER_DIR}/hiscanner_output/final_calls_lambda${lambda} ]; then
        echo "Running CNV calling with lambda=$lambda"
        #REPLACE lambda_value
        cat ${HISCANNER_DIR}/config_orig.yaml | grep -v "lambda_value" > config.yaml
        echo "lambda_value: ${lambda}" >> config.yaml
        conda run -n hiscanner_test hiscanner run --step cnv
        mv ${HISCANNER_DIR}/hiscanner_output/final_calls ${HISCANNER_DIR}/hiscanner_output/final_calls_lambda${lambda}
    fi
done

exit

mkdir -p logs_lambda




exit



# Initialize curr_job_id to an empty string; no dependency for the first iteration
curr_job_id=""

echo "CURRENT DIR: $(realpath .)"

for lambda_value in "${lambda_values[@]}"; do
    echo "Running with lambda_value=$lambda_value"

    #### First CNV run ####
    log_file1="logs_lambda/error_${lambda_value}_first.out"
    log_file1_err="logs_lambda/error_${lambda_value}_first.err"
    if [[ -z "$curr_job_id" ]]; then
        # First iteration: no dependency
        jobid1=$(sbatch --parsable \
        -o "$log_file1" \
        -e "$log_file1_err" \
        -p medium,park \
        -t 1-00:00:00 \
        -c 1 \
        --mem 6G \
        --wrap="
            # Update lambda_value in config.yaml
            sed -i \"s/^lambda_value: [0-9]\\+\\(\\.[0-9]\\+\\)\\?/lambda_value: $lambda_value/\" \"$config_file\"
            
            # Reset metadata.txt to original
            cp \"${metadata_file}.bak\" \"$metadata_file\"
            
            # Run the first CNV step
            module load gcc/14.2.0 conda/miniforge3/24.11.3-0 samtools
            conda run -n hiscanner_test hiscanner run --step cnv
        ")
        echo "Submitted first CNV job: $jobid1"
    else
        # Subsequent iterations: dependent on jobid2 from the previous iteration
        jobid1=$(sbatch --parsable \
        -o "$log_file1" \
        -e "$log_file1_err" \
        -p medium,park \
        -t 1-00:00:00 \
        -c 1 \
        --mem 6G \
        --dependency=afterany:$curr_job_id \
        --wrap="
            # Update lambda_value in config.yaml
            sed -i \"s/^lambda_value: [0-9]\\+\\(\\.[0-9]\\+\\)\\?/lambda_value: $lambda_value/\" \"$config_file\"
            
            # Reset metadata.txt to original
            cp \"${metadata_file}.bak\" \"$metadata_file\"
            
            # Run the first CNV step
            module load gcc/14.2.0 conda/miniforge3/24.11.3-0 samtools
            conda run -n hiscanner_test hiscanner run --step cnv
        ")
        echo "Submitted first CNV job (dependent on previous job): $jobid1"
    fi

    # Set curr_job_id to the current jobid1
    curr_job_id=$jobid1

    #### Second CNV run ####
    log_file2="logs_lambda/error_${lambda_value}_second.out"
    log_file2_err="logs_lambda/error_${lambda_value}_second.err"
    jobid2=$(sbatch --parsable \
        -o "$log_file2" \
        -e "$log_file2_err" \
        -p medium,park \
        -t 1-00:00:00 \
        -c 1 \
        --mem 6G \
        --dependency=afterany:$jobid1 \
        --wrap="
            # Remove bad cells from metadata after first job
            if [[ -f \"$log_file1\" ]]; then
                while read -r line; do
                    if [[ \"\$line\" == *\"Error processing cell\"* ]]; then
                        cell=\$(echo \"\$line\" | cut -d ':' -f1 | awk '{print \$NF}')
                        if [[ -n \"\$cell\" ]]; then
                            echo \"Removing cell \$cell from metadata\"
                            sed -i \"/\$cell/d\" \"$metadata_file\"
                        fi
                    fi
                done < \"$log_file1\"
            else
                echo \"Log file $log_file1 not found.\"
            fi

            # Clean up final_calls directory
            if [ -d 'hiscanner_output/final_calls' ]; then
                echo 'Cleaning up final_calls directory...'
                rm -r hiscanner_output/final_calls
            else
                echo 'Directory hiscanner_output/final_calls not found.'
            fi

            # Run second CNV analysis step
            module load gcc/14.2.0 conda/miniforge3/24.11.3-0 samtools
            conda run -n hiscanner_test hiscanner run --step cnv
            
            # Move results after the run
            if [ -d 'hiscanner_output/final_calls' ]; then
                mv hiscanner_output/final_calls \"hiscanner_output/final_calls_lambda$lambda_value\"
            else
                echo 'Directory hiscanner_output/final_calls not found. Skipping move.'
            fi
        ")

    echo "Submitted second CNV job: $jobid2"

    # Update curr_job_id to jobid2 for the next iteration
    curr_job_id=$jobid2

done

#mv hiscanner_output/final_calls hiscanner_output/final_calls_lambda0.25
    
#rm config.yaml metadata.txt
