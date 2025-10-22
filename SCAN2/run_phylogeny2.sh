#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH -p priopark
#SBATCH -A park_contrib

echo "making matrices for standard SCAN2 calls.."
standard=$(ls phylogeny/pass_set_102_cells_*_dnvs.txt)
outalt=phylogeny/scan2_alt_matrix_$(basename $standard .txt|cut -f5- -d_).csv
outdp=phylogeny/scan2_dp_matrix_$(basename $standard .txt|cut -f5- -d_).csv
if [ ! -f $outalt ] && [ ! -f ${outalt}.gz ] ; then
    phylogeny/make_sequoia_matrices.R $outalt $outdp phylogeny/subsetted_objects/*.txt
else
    echo "skipping $outalt"
fi
gzip $outalt $outdp


echo "making matrices for relaxed SCAN2 calls.."
relaxed=$(ls phylogeny/pass_set_102_cells_*_dnvs.relaxed.txt)
outalt=phylogeny/scan2_alt_matrix_$(basename $relaxed .txt|cut -f5- -d_).csv
outdp=phylogeny/scan2_dp_matrix_$(basename $relaxed .txt|cut -f5- -d_).csv
if [ ! -f $outalt ] && [ ! -f ${outalt}.gz ] ; then
    phylogeny/make_sequoia_matrices.R $outalt $outdp phylogeny/subsetted_objects_relaxed/*.txt
else
    echo "skipping $outalt"
fi
gzip $outalt $outdp

