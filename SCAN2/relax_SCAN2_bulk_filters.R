#!/usr/bin/env Rscript
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH -t 20:00
#SBATCH -p priopark
#SBATCH -A park_contrib

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop('usage: ./update_max_bulk_binom.R in.rda out.rda')
}

in.rda <- args[1]
out.rda <- args[2]

if (file.exists(out.rda))
    stop(paste('output file', out.rda, 'already exists, please delete it first'))

library(scan2)

load(in.rda, verb=TRUE)

# Allow for reads in bulk so long as the probability of being a germline mutation is very low.
relaxed.params=list(max.bulk.alt=100, max.bulk.af=1, max.bulk.binom.prob=1e-5)

results <- update.static.filter.params(results, new.params=list(snv=relaxed.params, indel=relaxed.params))

save(results, file=out.rda, compress=FALSE)
