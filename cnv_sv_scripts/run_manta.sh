#!/bin/bash

# Script to run Manta somatic structural variant calling
#
# This script configures and executes the Manta workflow using single cell PTA and normal
# BAM files along with a reference genome. The workflow is created inside the
# specified working directory and executed locally.
#
# Input:
#   - sc PTA BAM file
#   - Normal (blood) BAM file
#   - Reference FASTA file
#   - Working directory for the Manta run
#   - Optional: number of CPU cores to use (default: 8)
#
# Output:
#   - Manta workflow directory containing variant call results
#
# Dependencies:
#   - Manta

set -euo pipefail

if [ $# -lt 4 ] || [ $# -gt 5 ]; then
    echo "Usage: $0 <tumor_bam> <normal_bam> <reference_fasta> <working_dir> [cores]"
    echo "Example: $0 tumor.bam normal.bam reference.fa /path/to/workdir 12"
    exit 1
fi

sample_bam=$1
blood_bam=$2
reference_fasta=$3
working_dir=$4
cores=${5:-8}

# Ensure input files exist before proceeding
if [ ! -f "$sample_bam" ]; then
    echo "Error: Tumor BAM file '$sample_bam' not found"
    exit 1
fi

if [ ! -f "$blood_bam" ]; then
    echo "Error: Normal BAM file '$blood_bam' not found"
    exit 1
fi

if [ ! -f "$reference_fasta" ]; then
    echo "Error: Reference FASTA file '$reference_fasta' not found"
    exit 1
fi

# Create working directory and configure the workflow
mkdir -p "$working_dir"

configManta.py \
    --normalBam="$blood_bam" \
    --tumorBam="$sample_bam" \
    --referenceFasta="$reference_fasta" \
    --runDir="$working_dir"

"$working_dir/runWorkflow.py" -m local -j "$cores"
