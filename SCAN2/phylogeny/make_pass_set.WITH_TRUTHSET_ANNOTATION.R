#!/usr/bin/env Rscript
#SBATCH -c 1
#SBATCH --mem=12G
#SBATCH -t 15:00
#SBATCH -p priopark
#SBATCH -A park_contrib

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop('usage: ./subset_scan2_object.R in_muts.txt.gz out.txt')
}

in.txt <- args[1]
out.txt <- args[2]
#in.rdas <- args[-1]

if (file.exists(out.txt))
    stop(paste('output file', out.txt, 'already exists, please delete it first'))

suppressMessages(library(scan2))

# Lung truth set
ts.lung <- fread("MT_ST002-1D_Truthset_v1.txt")[, in.ts := 'lung']
ts.colon <- fread("MT_ST002-1G_Truthset_v1.txt")[, in.ts := 'colon']

ts <- merge(ts.lung, ts.colon,
    all.x=T, all.y=T,
    by=c('chr','pos','refnt','altnt','qual','filter','rgn','mosaic.type'),
    suffixes=c('.lung', '.colon'))
ts[, in.ts := ifelse(is.na(in.ts.lung), 'colon', ifelse(is.na(in.ts.colon), 'lung', 'both'))]
ts <- ts[, .(chr, pos, refnt, altnt, rgn, mosaic.type, in.ts, af.colon, lr.af.colon, af.lung, lr.af.lung)]

p <- fread(in.txt)

# Group across samples so one row = one site called in any sample
p <- p[,.(n.pass=.N, min.bulk.dp=min(bulk.dp)), by=.(chr,pos,refnt,altnt,muttype,mutsig)]
p <- ts[p, , on=.(chr,pos,refnt,altnt)] 
p[is.na(in.ts), in.ts := 'no']

fwrite(p, file=out.txt, sep='\t')
