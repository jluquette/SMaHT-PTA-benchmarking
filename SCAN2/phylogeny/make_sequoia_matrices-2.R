#!/usr/bin/env Rscript
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 15:00
#SBATCH -p priopark
#SBATCH -A park_contrib

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
    stop('usage: ./make_sequoia_matrices.R out_alt_matrix.txt out_dp_matrix.txt data1.txt [ data2.txt ... dataN.txt ]')
}

out.alt.txt <- args[1]
out.dp.txt <- args[2]
in.txts <- args[-(1:2)]

if (file.exists(out.alt.txt))
    stop(paste('output file', out.alt.txt, 'already exists, please delete it first'))
if (file.exists(out.dp.txt))
    stop(paste('output file', out.dp.txt, 'already exists, please delete it first'))

suppressMessages(library(data.table))

x <- rbindlist(lapply(in.txts, fread))

fwrite(dcast(x, chr+pos+dbsnp+refnt+altnt+muttype+mutsig ~ sample, value.var='scalt'), file=out.alt.txt)
fwrite(dcast(x, chr+pos+dbsnp+refnt+altnt+muttype+mutsig ~ sample, value.var='dp'), file=out.dp.txt)
