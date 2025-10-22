#!/bin/bash

# Script to perform CNV analysis using CNVpytor
# 
# This script processes BAM and phased VCF files to detect Copy Number Variations
# using CNVpytor tool. It generates both pytor and call files as output.
#
# Input:
#   - BAM file: Aligned sequence reads
#   - Phased VCF file: Phased variant calls
#   - Optional: number of CPU cores to use (default: 8)
#
# Output:
#   - .pytor file: CNVpytor binary file containing histogram and statistics
#   - .calls file: Contains CNV calls and their coordinates
#
# Dependencies:
#   - CNVpytor
#   - SAMtools

set -euo pipefail

# Check if correct number of arguments is provided
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <bam_file> <vcf_file> [cores]"
    echo "Example: $0 sample.bam sample.vcf 12"
    exit 1
fi

# Assign arguments to variables
bam=$1
vcf=$2
cores=${3:-8}

# Check if input files exist
if [ ! -f "$bam" ]; then
    echo "Error: BAM file '$bam' not found"
    exit 1
fi

if [ ! -f "$vcf" ]; then
    echo "Error: VCF file '$vcf' not found"
    exit 1
fi

# Extract sample name from BAM file (assuming format like sample.bam)
sample=$(basename "$bam" .bam)

# Run CNVpytor commands
cnvpytor -root ${sample}.pytor -log ${sample}.log -rd ${bam} -j ${cores} -chrom `seq -f "chr%g" 1 22` chrX chrY chrM
cnvpytor -root ${sample}.pytor -log ${sample}.log -snp ${vcf} -nofilter -chrom `seq -f "chr%g" 1 22` chrX chrY chrM
cnvpytor -root ${sample}.pytor -log ${sample}.log -mask_snps
cnvpytor -root ${sample}.pytor -log ${sample}.log -pileup ${bam} -j ${cores}
cnvpytor -root ${sample}.pytor -log ${sample}.log -his 10000 100000
cnvpytor -root ${sample}.pytor -log ${sample}.log -baf 10000 100000 -usephase
cnvpytor -root ${sample}.pytor -log ${sample}.log -call combined 10000 -usephase -mindbaf 0.2 > ${sample}.2d.10k.tsv
cnvpytor -root ${sample}.pytor -log ${sample}.log -call combined 100000 -usephase -mindbaf 0.2 > ${sample}.2d.100k.tsv
