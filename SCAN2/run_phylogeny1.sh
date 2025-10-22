#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH -t 1:00:00
#SBATCH -p priopark
#SBATCH -A park_contrib


outbase=phylogeny/pass_set_102_cells
dirbase=phylogeny/subsetted_objects

# Standard SCAN2 calls with usual bulk filters
tmpout=${outbase}_NMUTS_muts.txt

echo "making non-relaxed sites list.."
phylogeny/make_pass_set.WITH_TRUTHSET_ANNOTATION.R all_SCAN2_pass_snv_indel_dnv.csv.gz $tmpout

nsnv=$(awk 'BEGIN { FS="\t"; } $12 == "snv"' $tmpout | wc -l | cut -f1 -d\ )
nindel=$(awk 'BEGIN { FS="\t"; } $12 == "indel"' $tmpout | wc -l | cut -f1 -d\ )
ndnv=$(awk 'BEGIN { FS="\t"; } $12 == "dnv"' $tmpout | wc -l | cut -f1 -d\ )
sitesf=${outbase}_${nsnv}_snvs_${nindel}_indels_${ndnv}_dnvs.txt

mv -n $tmpout $sitesf
echo $sitesf

mkdir -p $dirbase
for f in objs/full/*.rda; do
    sm=$(basename $f .rda)
    sbatch phylogeny/subset_scan2_object.R $f $sitesf $dirbase/$sm.txt
done


# For bulk-relaxed results
echo "making relaxed sites list.."
tmpout=${outbase}_NMUTS_muts.relaxed.txt

phylogeny/make_pass_set.WITH_TRUTHSET_ANNOTATION.R all_SCAN2_pass_snv_indel_dnv.relaxed.csv.gz $tmpout

nsnv=$(awk 'BEGIN { FS="\t"; } $12 == "snv"' $tmpout | wc -l | cut -f1 -d\ )
nindel=$(awk 'BEGIN { FS="\t"; } $12 == "indel"' $tmpout | wc -l | cut -f1 -d\ )
ndnv=$(awk 'BEGIN { FS="\t"; } $12 == "dnv"' $tmpout | wc -l | cut -f1 -d\ )
sitesf=${outbase}_${nsnv}_snvs_${nindel}_indels_${ndnv}_dnvs.relaxed.txt

mv -n $tmpout $sitesf
echo $sitesf

mkdir -p ${dirbase}_relaxed
for f in objs/full_bulk_binom_prob_1e-5/*.rda; do
    sm=$(basename $f .rda)
    sbatch phylogeny/subset_scan2_object.R $f $sitesf ${dirbase}_relaxed/$sm.txt
done
