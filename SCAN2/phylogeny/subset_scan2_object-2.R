#!/usr/bin/env Rscript
#SBATCH -c 1
#SBATCH --mem=30G
#SBATCH -t 30:00
#SBATCH -p priopark
#SBATCH -A park_contrib

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
    stop('usage: ./subset_scan2_object.R in.rda sites.txt out.txt')
}

in.rda <- args[1]
sites.txt <- args[2]
out.txt <- args[3]

if (file.exists(out.txt))
    stop(paste('output file', out.txt, 'already exists, please delete it first'))

library(scan2)

sites <- fread(sites.txt)

load(in.rda) # loads 'results'

# this works and only takes ~30s or so, but the mapping isn't perfect for MNVs.
# in particular, if there are two adjacent non-ref bases in data(.), they will
# get merged into an MNV but if the first base is not PASSed in any single cell,
# then it will be treated as a SNV in passing(.) downstream of this script.
#
# although this isn't perfect, it is pretty close to capturing SCAN2's MNV calls:
#
#   > reference.passing.set.across.cells[, table(muttype)]
#   muttype
#      dbs  indel    mnv    snv 
#     1234   6766     49 114923 
#   > helper.combine.mnv(data(results))[, table(muttype)]
#   muttype
#      dbs  indel    mnv    snv 
#     1166   6766     42 114413
#
# About 95% of SCAN2's DBSes are recovered at a loss of 0.5% of SNVs.
ss <- helper.combine.mnv(data(results))
ss <- ss[muttype != 'mnv']  # Let's not worry about MNVs for now

# is.na(sample) shouldn't happen with properly analyzed data
ss <- ss[sites,,on=.(chr, pos, refnt, altnt, muttype, mutsig)][!is.na(sample)]
ss[in.ts == "", in.ts := "no"]


fwrite(ss, file=out.txt)
