#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --job-name=mpileup
#SBATCH --mem=2G
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=24
#SBATCH -p day


#You need to install BSMN-pipeline: https://github.com/bsmn/bsmn-pipeline as conda environment

module load miniconda
conda activate bp

path="/vast/palmer/pi/vaccarino/agn23/scWGS/collab"
input="varinats.txt" #file should have 4 columns: CHR POS REF ALT 
out="/vast/palmer/scratch/vaccarino/agn23/pileup"
mkdir -p ${out}
for sample in `cat ${path}/files.txt`;
do
python3 /gpfs/gibbs/project/vaccarino/agn23/bsmn_pipeline/utils/somatic_vaf.2.py -q 20 -Q 20 -b ${path}/${sample}.bam  -n 24 ${input} > ${out}/${sample}.txt
done

conda deactivate bp
module purge
module load R/4.3.0-foss-2020b
export R_LIBS=/gpfs/gibbs/project/vaccarino/agn23/R/4.3

R --vanilla <<code
library(tidyverse)
QC = read.csv("/gpfs/gibbs/pi/vaccarino/agn23/scWGS/collab/annotated_qc_fails.csv",header = TRUE)
fail = QC %>% filter(n.qc.fails > 1) %>% select(up_id_bam)
L = read.table("/vast/palmer/scratch/vaccarino/agn23/pileup/ST002-1D.txt",header = FALSE)
C = read.table("/vast/palmer/scratch/vaccarino/agn23/pileup/ST002-1G.txt",header = FALSE)
LV = L %>% mutate(join = paste0(V1,"_",V2,"_",V3,"_",V4)) %>% filter(V8 > 0) %>% select(join,V8) %>% rename(Bulk_Lung = V8)
CV = C %>% mutate(join = paste0(V1,"_",V2,"_",V3,"_",V4)) %>% filter(V8 > 0) %>% select(join,V8) %>% rename(Bulk_Colon = V8)
LR = L %>% mutate(join = paste0(V1,"_",V2,"_",V3,"_",V4)) %>% filter(V8 > 0) %>% select(join,V8) %>% rename(Bulk_Lung = V6)
CR = C %>% mutate(join = paste0(V1,"_",V2,"_",V3,"_",V4)) %>% filter(V8 > 0) %>% select(join,V8) %>% rename(Bulk_Colon = V6)
NR = merge(LR,CR,all = TRUE)
NV = merge(LV,CV,all = TRUE)
s = read.table(paste0("/gpfs/gibbs/pi/vaccarino/agn23/scWGS/collab/files.txt"), header = FALSE)
    for (sample in s$V1) {
    A = read.table(paste0("/vast/palmer/scratch/vaccarino/agn23/pileup/",sample,".txt"), header = FALSE)
    V = A %>% mutate(join = paste0(V1,"_",V2,"_",V3,"_",V4)) %>% filter(V8 > 0) %>% select(join,V8)
    R = A %>% mutate(join = paste0(V1,"_",V2,"_",V3,"_",V4)) %>% filter(V6 > 0) %>% select(join,V6)
    colnames(V)[2] = sample
    colnames(R)[2] = sample
    NR = merge(NR,R,all = TRUE)
    NV = merge(NV,V,all = TRUE)
    }
m = grep(paste(fail$up_id_bam, collapse = "|"),colnames(NR))
NR = NR[-m]
m = grep(paste(fail$up_id_bam, collapse = "|"),colnames(NV))
NV = NV[-m]
NR = data.frame(NR[4:ncol(NR)], row.names = NR$join)
NR %>% write.table("NR.txt",sep = "\t",quote = FALSE)
NV = data.frame(NV[4:ncol(NV)], row.names = NV$join)
NV %>% write.table("NV.txt",sep = "\t",quote = FALSE)
code

#files NV and NR should be used with Sequoia to create lineage tree