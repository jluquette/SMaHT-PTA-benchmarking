options(stringsAsFactors = F)
#----------------------------------
# Load packages (install if they are not installed yet)
#----------------------------------
options(stringsAsFactors = F)
cran_packages=c("ggplot2","ape","seqinr","stringr","data.table","tidyr","dplyr","VGAM","MASS","devtools")
bioconductor_packages=c("Rsamtools","GenomicRanges")

for(package in cran_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    install.packages(as.character(package),repos = "http://cran.us.r-project.org")
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if (!require("BiocManager", quietly = T, warn.conflicts = F))
  install.packages("BiocManager")
for(package in bioconductor_packages){
  if(!require(package, character.only=T,quietly = T, warn.conflicts = F)){
    BiocManager::install(as.character(package))
    library(package, character.only=T,quietly = T, warn.conflicts = F)
  }
}
if(!require("treemut", character.only=T,quietly = T, warn.conflicts = F)){
  install_git("https://github.com/NickWilliamsSanger/treemut")
  library("treemut",character.only=T,quietly = T, warn.conflicts = F)
}

#----------------------------------
# Functions
#----------------------------------

plot_spectrum = function(bed,save,add_to_title="",genomeFile = "/nfs/cancer_ref01/Homo_sapiens/37/genome.fa"){
  mutations=as.data.frame(bed)
  colnames(mutations) = c("chr","pos","ref","mut")
  mutations$pos=as.numeric(mutations$pos)
  mutations = mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")) & mutations$chr %in% c(paste0("chr",c(1:22,"X","Y")),c(1:22,"X","Y")),]
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(as.numeric(mutations$pos)-1, 
                                                                                     as.numeric(mutations$pos)+1))))
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs = table(paste(mutations$sub,paste(substr(mutations$trinuc_ref_py,1,1),substr(mutations$trinuc_ref_py,3,3),sep="-"),sep=","))
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_full = freqs[full_vec]; freqs_full[is.na(freqs_full)] = 0; names(freqs_full) = full_vec
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  if(!is.null(save)) pdf(save,width=12,height=4)
  if(is.null(save)) dev.new(width=12,height=4)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y = freqs_full; maxy = max(y)
  h = barplot(y, las=2, col=colvec, border=NA, ylim=c(0,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab="Number mutations", main=paste0("Number of mutations: ",sum(freqs_full), add_to_title))
  for (j in 1:length(sub_vec)) {
    xpos = h[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.2, xpos[2]+0.5, maxy*1.3, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.3, pos=3, label=sub_vec[j])
  }    
  if(!is.null(save)) dev.off()
  
}

exact.binomial=function(gender,NV,NR,cutoff=-5,qval_return=F){
  # Function to filter out germline variants based on unmatched
  # variant calls of multiple samples from same individual (aggregate coverage
  # ideally >150 or so, but will work with less). NV is matrix of reads supporting 
  # variants and NR the matrix with total depth (samples as columns, mutations rows, 
  # with rownames as chr_pos_ref_alt or equivalent). Returns a logical vector, 
  # TRUE if mutation is likely to be germline.
  
  XY_chromosomal = grepl("X|Y",rownames(NR))
  autosomal = !XY_chromosomal
  
  if(gender=="female"){
    NV_vec = rowSums(NV)
    NR_vec = rowSums(NR)
    pval = rep(1,length(NV_vec))
    for (n in 1:length(NV_vec)){
      if(NR_vec[n]>0){
        pval[n] = binom.test(x=NV_vec[n],
                             n=NR_vec[n],
                             p=0.5,alt='less')$p.value
      }
    }
  }
  # For male, split test in autosomal and XY chromosomal part
  if(gender=="male"){
    pval=rep(1,nrow(NV))
    NV_vec = rowSums(NV)[autosomal]
    NR_vec = rowSums(NR)[autosomal]
    pval_auto = rep(1,sum(autosomal))
    pval_XY = rep(1,sum(XY_chromosomal))
    
    for (n in 1:sum(autosomal)){
      if(NR_vec[n]>0){
        pval_auto[n] = binom.test(x=NV_vec[n],
                                  n=NR_vec[n],
                                  p=0.5,alt='less')$p.value
      }
    }
    NV_vec = rowSums(NV)[XY_chromosomal]
    NR_vec = rowSums(NR)[XY_chromosomal]
    for (n in 1:sum(XY_chromosomal)){
      if(NR_vec[n]>0){
        pval_XY[n] = binom.test(x=NV_vec[n],
                                n=NR_vec[n],
                                p=0.95,alt='less')$p.value
      }
    }
    pval[autosomal]=pval_auto
    pval[XY_chromosomal]=pval_XY
  }
  qval = p.adjust(pval,method="BH")
  if(qval_return){
    return(qval)
  }else{
    germline = log10(qval)>cutoff
    return(germline)
  }
}

estimateRho_gridml = function(NV_vec,NR_vec) {
  # Function to estimate maximum likelihood value of rho for beta-binomial
  rhovec = 10^seq(-6,-0.05,by=0.05) # rho will be bounded within 1e-6 and 0.89
  mu=sum(NV_vec)/sum(NR_vec)
  ll = sapply(rhovec, function(rhoj) sum(dbetabinom(x=NV_vec, size=NR_vec, rho=rhoj, prob=mu, log=T)))
  return(rhovec[ll==max(ll)][1])
}

beta.binom.filter = function(NR,NV){
  # Function to apply beta-binomial filter for artefacts. Works best on sets of
  # clonal samples (ideally >10 or so). As before, takes NV and NR as input. 
  # Optionally calculates pvalue of likelihood beta-binomial with estimated rho
  # fits better than binomial. This was supposed to protect against low-depth variants,
  # but use with caution. Returns logical vector with good variants = TRUE
  
  rho_est = pval = rep(NA,nrow(NR))
  for (k in 1:nrow(NR)){
    rho_est[k]=estimateRho_gridml(NV_vec = as.numeric(NV[k,]),
                                  NR_vec=as.numeric(NR[k,]))
  }
  return(rho_est)
}

dbinomtrunc = function(x, size, prob, minx=4) {
  dbinom(x, size, prob) / pbinom(minx-0.1, size, prob, lower.tail=F)
}

estep = function(x,size,p.vector,prop.vector,ncomp, mode){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    if(mode=="Truncated") p.mat_estep[,i]=prop.vector[i]*dbinomtrunc(x,size,prob=p.vector[i])
    if(mode=="Full") p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)   
}

em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust,binom_mode){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities 
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust,mode=binom_mode)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust,mode=binom_mode)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6, mode="Full"){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC) 
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol,binom_mode=mode)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}

binom_pval_matrix = function(NV,NR,gender,qval_return=F) {
  NR_nonzero=NR
  NR_nonzero[NR_nonzero==0]=1
  pval_mat <- matrix(0, nrow = nrow(NV), ncol = ncol(NV))
  rownames(pval_mat)=rownames(NV)
  colnames(pval_mat)=colnames(NV)
  if(gender == "male") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NV)) {
        if (!grepl("X|Y",rownames(NV)[1])) {pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.5, alternative = "less")$p.value}
        else {pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.95, alternative = "less")$p.value}
      }
    }
  } else if(gender == "female") {
    for(i in 1:nrow(NV)) {
      for (j in 1:ncol(NV)) {
        pval_mat[i,j] <- binom.test(NV[i,j], NR_nonzero[i,j], p = 0.5, alternative = "less")$p.value
      }
    }
  }
  if(qval_return){
    qval_mat=matrix(p.adjust(as.vector(pval_mat), method='BH'),ncol=ncol(pval_mat))
    rownames(qval_mat)=rownames(NV)
    colnames(qval_mat)=colnames(NV)
    return(qval_mat)
  }else{
    return(pval_mat)
  }
}


apply_mix_model=function(NV,NR,plot=T,min_clonal_mut_num=min_clonal_mut){
  peak_VAF=rep(0,ncol(NV))
  names(peak_VAF)=colnames(NV)
  autosomal=!grepl("X|Y",rownames(NV))
  for(s in colnames(NV)){
    muts_include=NV[,s]>3&autosomal
    if(sum(muts_include)>5){
      NV_vec=NV[muts_include,s]
      NR_vec=NR[muts_include,s]
      res=binom_mix(NV_vec,NR_vec,mode="Truncated",nrange=1:3)
      saveRDS(res,paste0(output_dir,s,"_binom_mix.Rdata"))
      
      if(plot){
        pdf(paste0(output_dir,s,"_binom_mix.pdf"))
        p=hist(NV_vec/NR_vec,breaks=20,xlim=c(0,1),col='gray',freq=F,xlab="Variant Allele Frequency",
               main=paste0(s,", (n=",length(NV_vec),")"))
        cols=c("red","blue","green","magenta","cyan")
        
        y_coord=max(p$density)-0.5
        y_intv=y_coord/5
        
        for (i in 1:res$n){
          depth=rpois(n=5000,lambda=median(NR_vec))
          sim_NV=unlist(lapply(depth,rbinom,n=1,prob=res$p[i]))
          sim_VAF=sim_NV/depth
          sim_VAF=sim_VAF[sim_NV>3]
          dens=density(sim_VAF)
          lines(x=dens$x,y=res$prop[i]*dens$y,lwd=2,lty='dashed',col=cols[i])
          y_coord=y_coord-y_intv/2
          text(y=y_coord,x=0.9,label=paste0("p1: ",round(res$p[i],digits=2)))
          segments(lwd=2,lty='dashed',col=cols[i],y0=y_coord+y_intv/4,x0=0.85,x1=0.95)
        }
        dev.off()
      }
      peak_VAF[s]=max(res$p[(res$prop*length(res$Which_cluster))>min_clonal_mut])
    }
  }
  return(peak_VAF)
}

add_ancestral_outgroup=function(tree,outgroup_name="Ancestral"){
  #This function adds the ancestral tip at the end
  tmp=tree$edge
  N=length(tree$tip.label)
  newroot=N+2
  renamedroot=N+3
  ancestral_tip=N+1
  tmp=ifelse(tmp>N,tmp+2,tmp)
  
  tree$edge=rbind(c(newroot,renamedroot),tmp,c(newroot,ancestral_tip))
  tree$edge.length=c(0,tree$edge.length,0)
  
  tree$tip.label=c(tree$tip.label,outgroup_name)
  tree$Nnode=tree$Nnode+1
  mode(tree$Nnode)="integer"
  mode(tree$edge)="integer"
  return(tree)
}

low_vaf_in_pos_samples = function(NR, NV, gender, define_pos = 3, qval_return = F) {
  pval=rep(0,nrow(NR))
  if(gender == "male") {
    for(n in 1:nrow(NR)) {
      NV_vec=NV[n,]
      NR_vec=NR[n,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos=NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos=NR_vec[which(NV_vec >= define_pos)]
        if (grepl("X|Y",rownames(NR)[n])) {
          pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.95, alt = "less")$p.value
        } else {
          pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
        }
      }
    }
  } else if(gender == "female") {
    for(n in 1:nrow(NR)) {
      NV_vec=NV[n,]
      NR_vec=NR[n,]
      if(any(NV_vec >= define_pos)){
        NV_vec_pos=NV_vec[which(NV_vec >= define_pos)]
        NR_vec_pos=NR_vec[which(NV_vec >= define_pos)]
        pval[n]=binom.test(sum(NV_vec_pos), sum(NR_vec_pos), p = 0.5, alt = "less")$p.value
      }
    }
  }
  if(qval_return){
    return(p.adjust(pval,method="BH"))
  }else{
    return(pval)
  }
}

#----------------------------------
# Read in data
#----------------------------------
NR=read.csv("~/Downloads/scan2_dp_matrix_132369_muts.relaxed.csv.gz")
NV=read.csv("~/Downloads/scan2_alt_matrix_132369_muts.relaxed.csv.gz")
meta_data=read.csv("~/Downloads/sample_metadata.csv")

bad_cells=c("SM002-Lu-4-310","SMACS3WFHOU9","SMACSDW7W9JM","ST002-1G-Colon-1-PTA","ST002-1G-Colon-3-PTA","ST002-1D-Lung-1-PTA","ST002-1D-Lung-3-PTA","ST002-1D-Lung-4-PTA",
            "CollabPTA_Yale_colon_11","SMACSDLFV1VU","SM002-Lu-12-873","SMACSBM69H6Q","SMACSIGNW5HE","SMACSH1WE7O7")

Mut_ID=paste0(NV$chr,"_",NV$pos,"_",NV$ref,"_",NV$alt)
rownames(NV)=rownames(NR)=Mut_ID

NR=NR[,7:ncol(NR)]
NV=NV[,7:ncol(NV)]
colnames(NR)=colnames(NV)=gsub("\\.","-",colnames(NR))
NR=NR[,!colnames(NR)%in%bad_cells]
NV=NV[,!colnames(NV)%in%bad_cells]

NR_nonzero=NR
NR_nonzero[NR_nonzero==0]=1
VAF=NV/NR_nonzero

#Depth_filter = rowMeans(NR,na.rm=T)>10&rowMeans(NR,na.rm=T)<100
#----------------------------------
# Filter mutations
#----------------------------------

XY_chromosomal = grepl("X|Y",Mut_ID)
autosomal = !XY_chromosomal
xy_depth=mean(rowMeans(NR[XY_chromosomal,]),na.rm=T)
autosomal_depth=mean(rowMeans(NR[autosomal,]),na.rm=T)

gender='male'
min_cov=10
max_cov=100
if(xy_depth>0.8*autosomal_depth) gender='female'
Depth_filter = (rowMeans(NR,na.rm = T)>min_cov&rowMeans(NR,na.rm = T)<max_cov&autosomal)|
  (rowMeans(NR,na.rm = T)>(min_cov/2)&rowMeans(NR,na.rm = T)<(max_cov/2)&XY_chromosomal)


centers=meta_data$center[meta_data$amp!="bulk"]
names(centers)=meta_data$sample[meta_data$amp!="bulk"]
centers=centers[gsub("\\.","-",colnames(NR))]
rho_table=data.frame(Mut=Mut_ID)
rownames(rho_table)=rho_table$Mut
unique_centers=unique(centers)
for(center in unique_centers){
  center_samples=centers==center
  rho_est_shared=beta.binom.filter(NR=NR[,center_samples], NV=NV[,center_samples])
  rho_table[,paste0("Rho_",center)]=rho_est_shared
  rho_table[,paste0("Present_",center)]=rowSums(NV[,center_samples]>2)|rowSums(VAF[,center_samples]>0.15)
  print(center)
}
rho_table$rho_overall=beta.binom.filter(NR=NR, NV=NV)

rho_table$Over_0.1_present=rho_table$Over_0.35_present=0
for(n in 1:nrow(rho_table)){
  rho_table$Over_0.1_present[n]=sum(rho_table[n,grepl("Rho",colnames(rho_table))]>0.1&rho_table[n,grepl("Present",colnames(rho_table))]>0)
  rho_table$Over_0.35_present[n]=sum(rho_table[n,grepl("Rho",colnames(rho_table))]>0.35&rho_table[n,grepl("Present",colnames(rho_table))]>0)
}
rho_table$Total_present=rowSums(rho_table[,grepl("Present",colnames(rho_table))]>0)

for(n in 1:nrow(rho_table)){
  rho_table$Good[n]=rowSums(NV[n,]>2)==1&rowSums(VAF[n,]>0.15)==1
  if(!rho_table$Good[n]&!is.na(rho_table$Good[n])){
    if(rho_table$Total_present[n]==1&!is.na(rho_table$Total_present[n])){
      rho_table$Good[n]=rho_table$rho_overall[n]>0.1&rho_table$Over_0.2_present[n]>0
    }else{
      threshold=rho_table$Total_present[n]
      rho_table$Good[n]=rho_table$Over_0.35_present[n]>=threshold
    } 
  }
}

bad_muts=Mut_ID[Mut_ID%in%c(rho_table$Mut[!rho_table$Good|is.na(rho_table$Good)],names(Depth_filter)[!Depth_filter])]

NV_filtered=NV[!rownames(NV)%in%bad_muts,]
NR_filtered=NR[!rownames(NR)%in%bad_muts,]
write.table(NV_filtered,"~/Desktop/SMaHT_PTA/July/NV_filtered.txt")
write.table(NR_filtered,"~/Desktop/SMaHT_PTA/July/NR_filtered.txt")

NR_flt_nonzero=NR_filtered
NR_flt_nonzero[NR_flt_nonzero==0|is.na(NR_flt_nonzero)]=1
NV_filtered_nonNa=NV_filtered
NV_filtered_nonNa[is.na(NV_filtered_nonNa)]=0

#----------------------------------
# Make genotype matrix for tree
#----------------------------------

genotype_bin=as.matrix(NV_filtered_nonNa/NR_flt_nonzero)

genotype_bin[genotype_bin<0.1]=0
genotype_bin[genotype_bin>=0.25]=1
genotype_bin[genotype_bin>0&genotype_bin<1]=0.5
present_vars_full=rowSums(genotype_bin>0)>0
colnames(NV_filtered)=colnames(NR_filtered)=colnames(genotype_bin)=gsub("scalt_","",colnames(genotype_bin))

Ref = rep("A",nrow(genotype_bin))
Alt = rep("T",nrow(genotype_bin))
dna_strings = list()
dna_strings[1]=paste(Ref,sep="",collapse="") #Ancestral sample
for (n in 1:ncol(genotype_bin)){
  Mutations = Ref
  Mutations[genotype_bin[,n]==0.5] = '?'
  Mutations[genotype_bin[,n]==1] = Alt[genotype_bin[,n]==1]
  dna_string = paste(Mutations,sep="",collapse="")
  dna_strings[n+1]=dna_string
}

names(dna_strings)=c("Ancestral",colnames(genotype_bin))
require(seqinr)
output_dir="~/Desktop/SMaHT_PTA/July/"
patient_ID="Benchmark_v2"
mut_id="both"
write.fasta(dna_strings, names=names(dna_strings),paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa"))

# run mpboot to get max pars tree
# system(paste0(path_to_mpboot,"mpboot -s ",output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa -bb 1000"),ignore.stdout = T)
#

#----------------------------------
# Map mutations to tree
#----------------------------------
tree=read.tree(paste0(output_dir,patient_ID,"_",mut_id,"_for_MPBoot.fa.treefile"))
tree=drop.tip(tree,"Ancestral")
tree$edge.length = rep(1, nrow(tree$edge)) 

NR_tree=NR_filtered[present_vars_full,]
NR_tree[is.na(NR_tree)]=30
NR_tree$Ancestral=30
NV_tree=NV_filtered[present_vars_full,]
NV_tree[is.na(NV_tree)]=0
NV_tree$Ancestral=0

p.error = rep(0.3,ncol(NV_tree))
p.error[colnames(NV_tree)=="Ancestral"]=1e-6
res=assign_to_tree(tree,
                   mtr=as.matrix(NV_tree),
                   dep=as.matrix(NR_tree),
                   error_rate = p.error)
tree_mut_pval=0.01
badvars=res$summary$p_else_where>tree_mut_pval
vaf=NV_tree/NR_tree
heatmap(t(vaf[badvars&rowSums(vaf>0.1)>1,]),scale='none')

heatmap(t(vaf[early_muts,]),scale='none')


NR_tree_rescue=NR[!rownames(NR)%in%rownames(NR_tree),]
NR_tree_rescue[is.na(NR_tree_rescue)]=30
NR_tree_rescue$Ancestral=30
NV_tree_rescue=NV[!rownames(NV)%in%rownames(NV_tree),]
NV_tree_rescue[is.na(NV_tree_rescue)]=0
NV_tree_rescue$Ancestral=0

p.error = rep(0.3,ncol(NV_tree_rescue))
p.error[colnames(NV_tree_rescue)=="Ancestral"]=1e-6
res_rescue=assign_to_tree(tree,
                          mtr=as.matrix(NV_tree_rescue),
                          dep=as.matrix(NR_tree_rescue),
                          error_rate = p.error)

tree_mut_pval=0.001
goodvars=res_rescue$summary$p_else_where<tree_mut_pval
vaf=NV_tree_rescue/NR_tree_rescue

heatmap(as.matrix(vaf[sample(x=which(goodvars),size=2000),]),scale='none')


png("~/Desktop/SMaHT_PTA/earliest_muts_heatmap.png", width = 1000, height = 1000)
heatmap(as.matrix(vaf[earliest_muts,]),scale='none')
dev.off()

edge_length_nonzero = table(res$summary$edge_ml[res$summary$p_else_where<tree_mut_pval])
edge_length = rep(0,nrow(tree$edge))
names(edge_length)=1:nrow(tree$edge)
edge_length[names(edge_length_nonzero)]=edge_length_nonzero
tree$edge.length=as.numeric(edge_length)

saveRDS(res,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_tree_v2.Rdata"))
write.tree(tree, paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length_v2.tree"))

pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_branch_length_final.pdf"),height=15,width=7)
plot(tree)
axisPhylo(side = 1,backward=F)
dev.off()


#----------------------------------
# Final plots
#----------------------------------

tree_collapsed=di2multi(tree)
tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
tree_collapsed=drop.tip(tree_collapsed,"Ancestral")
pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_equal_branch_length_final.pdf"),height=15,width=7)
plot(tree_collapsed)
dev.off()

tree_collapsed_df=as.data.frame(fortify(tree_collapsed))
pdf(paste0(output_dir,patient_ID,"_",mut_id,"_tree_with_equal_branch_length.pdf"),height=15,width=7)

tree_df$tissue="Co"
tree_df$tissue[tree_df$label%in%meta_data$sample[meta_data$tissue=="lung"]]="Lu"

pdf("~/Desktop/SMaHT_PTA/July/tree_with_equal_branch_length_coloured.pdf",height=11,width=4)
plot(tree_collapsed,show.tip.label=F)
for(n in 1:sum(tree_collapsed_df$isTip)){
  if(grepl("Co",tree_df$tissue[n])){
    shape=21
  }else{
    shape=22
  }
  points(x=tree_collapsed_df$x[n],y=node.height(tree)[n],pch=shape,bg=center_cols[centers[tree_collapsed_df$label[n]]])
}
dev.off()

all_cols=c("grey80","firebrick","peachpuff","steelblue","salmon","gold","olivedrab3")
names(all_cols)=sigs
pdf("~/Desktop/SMaHT_PTA/tree_with_branch_length_coloured_v2.pdf",height=11,width=4)
plot(tree,show.tip.label=F)
# for(n in 1:sum(tree_collapsed_df$isTip)){
#   if(grepl("Co",tree_df$label[n])){
#     shape=21
#   }else{
#     shape=22
#   }
#   points(x=tree_df$x[n],y=node.height(tree)[n],pch=shape,bg=center_cols[centers[tree_df$label[n]]])
# }
axisPhylo(side = 1,backward=F)
#dev.off()


n=as.numeric(branches[k])
x_end=tree_df$x[n]
x_start=tree_df$x[tree_df$parent[n]]
x_intv=x_end-x_start
y=node.height(tree)[n]
tipnum=sum(tree_df$isTip)
for (s in sigs){
  x_end=x_start+exposures[s,samples[k]]*x_intv
  rect(ybottom=y-min(0.015*tipnum,0.3),ytop=y+min(0.015*tipnum,0.3),xleft=x_start,xright=x_end,col=all_cols[s],lwd=0.25)
  x_start=x_end
}
dev.off()

Mutations_per_branch=as.data.frame(matrix(ncol=4,unlist(strsplit(rownames(NR_tree),split="_")),byrow = T))
rownames(Mutations_per_branch)=rownames(NR_tree)
colnames(Mutations_per_branch)=c("Chr","Pos","Ref","Alt")
Mutations_per_branch$Branch = tree$edge[res$summary$edge_ml,2]
Mutations_per_branch=Mutations_per_branch[res$summary$p_else_where<tree_mut_pval,]
Mutations_per_branch$Patient = patient_ID
Mutations_per_branch$SampleID = paste(patient_ID,Mutations_per_branch$Branch,sep="_")
write.table(Mutations_per_branch,paste0(output_dir,patient_ID,"_",mut_id,"_assigned_to_branches.txt"),quote=F,row.names=F,sep="\t")


rownames(Mutations_per_branch)=paste0(Mutations_per_branch$Chr,"_",Mutations_per_branch$Pos,"_",Mutations_per_branch$Ref,"_",Mutations_per_branch$Alt)

find_children = function(node, tree=tree_df){
  child_nodes = tree$node[tree$parent==node]
  tips = c()
  for (k in 1:length(child_nodes)){
    if (tree$isTip[tree$node==child_nodes[k]]){
      tips=c(tips,child_nodes[k])
    }
    else{
      tips=c(tips,find_children(child_nodes[k]))
    }
  }
  return(tips)
}

tree_df=as.data.frame(fortify(tree))
tree_df$samples_names=tree_df$label
for(n in which(!tree_df$isTip&tree_df$node!=tree_df$parent)){
  tree_df$samples_names[n]=paste(tree_df$label[find_children(node=n,tree=tree_df)],collapse="; ")
}
Mutations_per_branch$SampleID_full=tree_df[Mutations_per_branch$Branch,"samples_names"]

early_muts=rownames(Mutations_per_branch)[Mutations_per_branch$Branch%in%which(!tree_df$isTip)]
earliest_muts=early_muts[rowSums(vaf[early_muts,]>0.15,na.rm = T)>3]
Mutations_per_branch$SampleID_full=tree_df[Mutations_per_branch$Branch,"samples_names"]

pdf("~/Desktop/SMaHT_PTA/July/short_tree.pdf",height=11,width=4)
plot(tree,x.lim=c(0,75))
axisPhylo(side = 1,backward=F,at=seq(0,100,by=10))

dev.off()

