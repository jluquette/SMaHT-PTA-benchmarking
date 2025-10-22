# Figure 2

# Exclude QC-fail cells in all panels of this figure
pcs <- meta[per.cell.stats,,on=.(sample)][n.fails<2]

# Read in duplex data
#dm <- fread("Duplex_group_Metrics_All_GCC_Tissue_Homogenate.xlsx_All_together.csv")
dm <- fread("Duplex_group_Metrics_All_GCC_Tissue_Homogenate.09162025.csv")
dmdt <- melt(as.data.table(t(dm[c(3:5,33,39),c(19:24, 26:31)]))[, c('tech', 'Tissue', 'snv_rate.per.bp', 'indel_rate.per.bp') := list(paste(V1, V2), gsub(pattern='ST002-1[DG]\n ', '', V3), as.numeric(V4), as.numeric(V5))][, .(tech, Tissue, snv_rate.per.bp, indel_rate.per.bp)], id.vars=1:2, measure.vars=measure(muttype, dummy, sep='_'), value.name='rate.per.bp')
cat("WARNING: excluding NanoSeq Hpy indel burdens due to unresolved low calls\n")
dmdt <- dmdt[!(tech == 'NanoSeq Hpy' & muttype == 'indel')]

# % genome analyzble per mut type
pb <- ggplot(pcs[variable=='fraction.genome.analyzable' & muttype %in% c('snv', 'indel')],
    aes(x=fct_shift(fct_recode(muttype, SNV='snv', Indel='indel')), y=100*value)) +
    geom_boxplot(outliers=FALSE, linewidth=0.25) +
    geom_jitter(width=0.2, size=1/4) +
    base.theme + manuscriptfont +
    xlab(element_blank()) + ylab("Percent genome callable")

# BURDENS -------------------- SNV and indel only (DNV burden not calculated)
# below is in units of muts per genome (autosomes only)
#z <- meta[melt(as.data.table(mutburden(sm), keep.rownames='sample'), id.vars=1, variable.name='muttype'),,on=.(sample)][order(value)]
# burdens as rates per bp
z <- meta[melt(as.data.table(mutburden(sm)/get.gbp.by.genome(sm[[1]])/1e9, keep.rownames='sample'), id.vars=1, variable.name='muttype'),,on=.(sample)][order(value)]
for (mt in c('snv', 'indel'))
    for (tt in c('lung', 'colon'))
        z[n.fails<2 & tissue == tt & muttype==mt, plot_x := .I]
dmdt <- dmdt[z[n.fails<2, .(width=.N, x=mean(plot_x)), by=.(Tissue, muttype)],,on=.(muttype,Tissue)]
cat('Comparison of PTA mut rates and duplex:\n')
print(z[n.fails<2, .(pta.median=median(value)), by=.(Tissue,muttype)][dmdt[, .(duplex.median=median(rate.per.bp), duplex.min=min(rate.per.bp), duplex.max=max(rate.per.bp)), by=.(Tissue,muttype)],,on=.(Tissue,muttype)][, median.reldiff := (pta.median-duplex.median)/duplex.median][, max.reldiff := (pta.median-duplex.max)/duplex.max][, duplex.range := (duplex.max-duplex.min)/duplex.min][])

pc <- ggplot(z[n.fails<2], aes(x=plot_x, y=value)) +
    stat_summary(data=dmdt, geom='crossbar', aes(width=width, x=x, y=rate.per.bp), inherit.aes=FALSE, fun=median, fun.max=max, fun.min=min, fill='lightblue', col=NA) +
    geom_hline(data=dmdt[, median(rate.per.bp), by=.(Tissue,muttype)], aes(yintercept=V1), linewidth=3/4, col='lightblue4') +
    geom_hline(data=z[n.fails<2, median(value), by=.(Tissue,muttype)], aes(yintercept=V1), linewidth=3/4, col='firebrick2') +
    geom_point(size=1/4) +
    scale_y_log10(labels=label_log()) + 
    scale_x_continuous(expand=expansion(c(0,0))) +
    ylab("Somatic mutation rate per cell (muts/bp)") +
    #facet_grid(~interaction(fct_recode(tissue, Colon='colon', Lung='lung'), fct_recode(factor(muttype, levels=c('snv', 'indel', 'dnv', 'mnv')), SNV='snv', Indel='indel', DNV='dnv', MNV='mnv'), sep=' '), scale='free_x') +
    facet_nested(~fct_recode(factor(muttype, levels=c('snv', 'indel')), SNV='snv', Indel='indel') + Tissue, scale='free_x', nest_line=element_line(color="black")) +
    xlab("Single cell") +
    base.theme + manuscriptfont +
    theme(strip.text=element_text(margin=margin(1/8,1/8,1/8,1/8, unit='line')),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(1/2, 'lines'),
        plot.margin=margin(b=0))


# CALLS ------------------------ SNV, indel, DNV
pd <- ggplot(meta[p[muttype != 'mnv',.(value=.N),by=.(sample,muttype)],,on=.(sample)][n.fails<2], aes(x=Tissue, y=value)) +
        geom_boxplot(outliers=FALSE, linewidth=0.25) +
        geom_jitter(width=0.2, size=1/4) +
        facet_wrap(~fct_recode(factor(muttype, levels=c('snv', 'indel', 'dnv', 'mnv')), SNV='snv', Indel='indel', DNV='dnv', MNV='mnv'), nrow=1, scale='free_y') +
        xlab(element_blank()) + ylab('Mutation calls per cell') +
        base.theme + manuscriptfont + theme(strip.text=element_text(margin=margin(0,0,0,0)), panel.spacing=unit(1/5, 'lines'), plot.margin=margin(0,0,0,0))

fp.tab <- dcast(pcs[variable %in% c('fps', 'calls')], Tissue + muttype ~ variable, value.var='value', fun.aggregate=sum)
fp.tab[muttype %in% c('dnv', 'mnv'), fps := 0]
fp.tab[, catalog.fdr := fps/calls]
fp.tab[, tps := calls - fps]

pe <- ggplot(melt(fp.tab[muttype != 'mnv',.(Tissue, muttype, tps, fps)], id.vars=c(1:2)), aes(x=Tissue)) +
        geom_bar(aes(weight=value), fill='black') +
        #geom_bar(aes(weight=value, fill=variable)) +
        #geom_label(data=fp.tab[muttype %in% c('snv', 'indel')], aes(x=tissue, y=fps, label=paste0(round(100*catalog.fdr), '%')), nudge_y=100, col='#C04657', size=7, size.unit='pt') +
        xlab(element_blank()) + ylab("Total mutation calls") +
        facet_wrap(~fct_recode(factor(muttype, levels=c('snv', 'indel', 'dnv', 'mnv')), SNV='snv', Indel='indel', DNV='dnv', MNV='mnv'), nrow=1, scale='free_y') +
        base.theme + manuscriptfont + theme(plot.margin=margin(0,0,0,0)) +
        theme(panel.spacing=unit(1/3, 'lines'), strip.text=element_text(margin=margin(0,0,0,0))) +
        #scale_fill_manual(values=c(tps='black', fps='#C04657')) +
        guides(fill='none') #+ labs(title=mt)

# Shared call numbers
zsh <- pr[,.N,by=.(chr,pos,refnt,altnt,muttype)][N>1]  
pg <- ggplot(zsh[muttype!='mnv'], aes(x=N)) +
    geom_line(stat='count', linewidth=0.2) +
    geom_point(stat='count', size=1/4) +
    facet_grid(fct_recode(factor(muttype, levels=c('snv', 'indel', 'dnv', 'mnv')), SNV='snv', Indel='indel', DNV='dnv', MNV='mnv')~., scale='free_y') +
    base.theme + manuscriptfont + theme(panel.spacing=unit(1/3, 'lines'), strip.text=element_text(margin=margin(0,0,0,0))) +
    ylab('Number of somatic mutations') + xlab("Number of cells sharing the mutation")

# Mutation spectra - normalized so each sample contributes a total weight of 1
zms <- melt(meta[,.(Tissue,sample)][p[muttype=='snv', as.list(table(sbs96(mutsig))), by=.(sample)],,on=.(sample)], id.vars=1:2)[, .(mutsig=variable, value=value/sum(value)) , by=.(Tissue,sample)]
ph <- ggplot(zms, aes(x=mutsig, weight=value, fill=mutsig)) +
    geom_sbs96() +
    facet_wrap(vars(Tissue), scale='free_y', ncol=1) +
    base.theme + manuscriptfont +
    ylab(element_blank()) + xlab(element_blank()) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.spacing=unit(1/2, units='line')) #, strip.text=element_text(margin=margin(0,0,0,0)))


# Duplex sequencing spectra
# XXX: Raw - are these correct/current?
# Aug. 12 2025 - no, they weren't
#tn <- fread('duplex_raw_trinucleotide_matrix.cleaned.txt')
# This is current as of 8/12/25:
tn <- fread('raw_trinucleotide_matrix_2025_07_25.txt')  # warns column 1 header is blank, which is true
colnames(tn)[1] <- 'sample'
dn <- melt(tn[grep(pattern='ST002_1[DG]', sample), c(sample='Duplex', lapply(.SD, sum)), .SDcols=-1, by=.(Tissue=ifelse(grepl(pattern='ST002_1D', sample), 'Lung', 'Colon'))], id.vars=1:2, variable.name='mutsig')[, value := as.numeric(value)][, value := value/sum(value), by=Tissue][]
dn[, mutsig := paste0(substr(mutsig,5,5), substr(mutsig,1,1), substr(mutsig,7,7), ':', substr(mutsig,1,3))]

# undoing the each sample sums to 1 normalization above
# does end up looking better - especially the C>A in lung
zms <- melt(meta[,.(n.fails,Tissue,sample)][p[muttype=='snv', as.list(table(sbs96(mutsig))), by=.(sample)],,on=.(sample)][n.fails<2][, n.fails:=NULL], id.vars=1:2, variable.name='mutsig')
ph <- ggplot(rbind(zms,dn), aes(x=mutsig, weight=value, fill=mutsig)) +
    geom_sbs96() +
    facet_grid(sample ~ Tissue, scale='free') +
    base.theme + manuscriptfont +
    ylab(element_blank()) + xlab(element_blank()) +
    theme(aspect.ratio=1/3, axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.spacing=unit(1/2, units='line')) #, strip.text=element_text(margin=margin(0,0,0,0)))

# All sigs in one column to make space for heatmap
zms2 <- zms[,sum(value), by=.(Tissue, mutsig)][, .(mutsig, sample='Single-cell', value=V1/sum(V1)), by=Tissue]
#ph_col <- ggplot(rbind(zms2, dn), aes(x=mutsig, weight=round(100*value,1), fill=mutsig)) +
# *in one row
ph_row <- ggplot(rbind(zms2, dn), aes(x=mutsig, weight=value, fill=mutsig)) +
    geom_sbs96() +
    facet_nested(~Tissue + sample, scale='free') +
    base.theme + manuscriptfont +
    ylab('Fraction of mutations') + xlab(element_blank()) +
    theme(aspect.ratio=1/3, axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    #theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing=unit(1/2, units='line'),
        ggh4x.facet.nestline=element_line(color="black"),
        strip.text=element_text(margin=margin(1/8,1/8,1/8,1/8, unit='line')))
# for printing
tomat <- function(x) {tmp <- dcast(x, Tissue+sample~mutsig); ret <- t(as.matrix(tmp[,-(1:2)])); colnames(ret) <- tmp[,paste(Tissue,sample)]; ret }
cos.sim.fun <- function(a) crossprod(a)/tcrossprod(sqrt(colSums(a^2)))
cat('Cosine similarity between PTA and duplex:\n')
print(cos.sim.fun(tomat(rbind(dn, zms[,.(value=sum(value)), by=.(Tissue, mutsig)][, sample := 'single cell']))))


# Compute per-cell cosine similarity to duplex mean
pf2i <- ggplot(zms[dn[, .(Tissue,mutsig, duplex=value)],,on=.(Tissue, mutsig)][, .(cos.sim=sum(value*duplex)/(sqrt(sum(value^2))*sqrt(sum(duplex^2)))),by=.(Tissue,sample)],
    aes(x=Tissue, y=cos.sim)) +
    geom_boxplot(outliers=FALSE, linewidth=0.25) +
    geom_jitter(width=0.2, size=1/4) +
    base.theme + manuscriptfont +
    xlab(element_blank()) + ylab('Cosine sim. per cell')

########## Example of a signature plot where bars represent the median signature ##########
########## and points are added per-channel to demonstrate the variability of    ##########
########## individual cells.                                                     ##########
# XXX: Unfortunately the APOBEC peaks are so high that it is hard to see the main
# signature, making the plot hard to compare to the duplex average.
#ggplot(zms[, median(value), by=mutsig], aes(x=mutsig, weight=V1, fill=mutsig)) +
#    geom_sbs96() +
#    geom_jitter(data=zms, aes(x=mutsig, weight=1, y=value), size=1/2, width=0.2) +
#    facet_grid(tissue~., scale='free')



### Correlation heatmap on SNV spectra
zhm <- dcast(zms, Tissue+sample ~ mutsig)
# zmat <- as.matrix(zhm[,-1], rownames=1)
# zmeta <- meta[zhm[,.(sample)],.(tissue,center,alias),on=.(sample)]
# setDF(zmeta, rownames=zmeta[,paste(center, alias)])
# cos.sim <- tcrossprod(zmat) / tcrossprod(sqrt(rowSums(zmat^2)))
# zmeta2 <- meta[zhm[sample != 'ST002-1G-Colon-4-PTA',.(sample)],.(tissue,center,alias),on=.(sample)]
# setDF(zmeta2, rownames=zmeta2[,paste(center, alias)])
# zmat2 <- as.matrix(zhm[sample != 'ST002-1G-Colon-4-PTA'][,-(1:2)], rownames=rownames(zmeta2))
# cos.sim2 <- tcrossprod(zmat2) / tcrossprod(sqrt(rowSums(zmat2^2)))
#
# 3: row names with centers are too long for figure. must use new short alias
zmeta3 <- meta[zhm[sample != 'ST002-1G-Colon-4-PTA',.(sample)],,on=.(sample)]
#zmeta3 <- meta[zhm[sample != 'SMACSWQEUWU3',.(sample)],,on=.(sample)]
setDF(zmeta3, rownames=zmeta3[,alias2])
zmat3 <- as.matrix(zhm[sample != 'ST002-1G-Colon-4-PTA'][,-(1:2)], rownames=rownames(zmeta3))
#zmat3 <- as.matrix(zhm[sample != 'SMACSWQEUWU3'][,-(1:2)], rownames=rownames(zmeta3))
#zhm4 <- zhm[ n.fails<2]
#ggcorrhm((apply(as.matrix(zhm4[,-(1:12)], rownames=zhm4[,paste(center,alias)]), 1, function(row) row/sum(row))), annot_rows_df=as.data.frame(zhm4[,c(1,10,12)], row.names=zhm4[,paste(center,alias)]), col_scale='D', annot_rows_params=list(size=3), names_x_side='bottom', cluster_rows=T, cluster_cols=T, limits=c(0.1,1), layout='bottomleft', include_diag=F, dend_height=3)

zzz <- ggcorrhm(t(zmat3), annot_rows_df=zmeta3[,c('Tissue', 'Dataset'),drop=F],
#zzz <- ggcorrhm(t(zmat4), annot_rows_df=zmeta4[,1:2,drop=F],
    col_scale='D',
    annot_rows_params=list(size=3),
    #annot_rows_col=list(Tissue=scale_fill_manual(values=tissue.colors), Dataset=scale_fill_manual(values=center.colors)),
    # temporary: hide legends
    annot_rows_col=list(Tissue=scale_fill_manual(values=tissue.colors, guide='none'),
        Dataset=scale_fill_manual(values=center.colors, guide='none')),
    cluster_rows=T, cluster_cols=T,
    limits=c(min(cor(t(zmat3))),1),
    layout='topleft',
    include_diag=T,
    #dend_height=1,
    dend_height=4,
    show_names_diag=F,
    #show_names_y=T,
    show_names_x=F,
    show_names_y=F,
    #show_dend_rows=FALSE,
    show_dend_cols=FALSE,
    return_data=T,
    cluster_method='complete', # 'ward.D2',
    annot_rows_names_side='bottom',
    #annot_cols_names_side='left',
    annot_names_size=2,
# use these to visually confirm the clusters
#split_rows=6,
#split_cols=6,
    legend_order=NA, # temporary: hide the legends
    border_col=NA)
pf2_heatmap <- zzz$plot + base.theme + manuscriptfont +
    xlab(element_blank()) + ylab(element_blank()) +
    theme(#legend.position='bottom',
        panel.spacing=unit(1/2, unit='mm'),  # for split_{rows,cols}
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank() ) #,
        #plot.margin=margin(0,0,0,0))  #, panel.spacing=unit(1, unit='mm'))
pf2_heatmap

# plot signatures for the top 5 clusters
# extract center/alias -> cluster mapping
# below is for zmeta2
#clmap <- meta[setnames(as.data.table(cutree(zzz$col_clustering, k=5), keep.rownames=T), c('rn', 'cluster'))[, c('center', 'alias') := list(sapply(strsplit(rn, ' '), `[`, 1), sapply(strsplit(rn, ' '), `[`, 2))],,on=alias2]
clmap <- meta[setnames(as.data.table(cutree(zzz$col_clustering, k=6), keep.rownames=T), c('alias2', 'cluster')),,on=.(alias2)]

# match the order of the signatures to the heatmap
cluster.order <- paste0('Cluster ', rev(unique(cutree(zzz$col_clustering,k=6)[zzz$col_clustering$order])))

#pf2_hmsigplot <- ggplot(p[clmap,,on=.(sample)][muttype=='snv'][,mutsig := sbs96(mutsig)],
pf2_hmsigplot <- ggplot(p[clmap,,on=.(sample)][muttype=='snv'][,mutsig := sbs96(mutsig)][,.(n=.N),by=.(cluster,sample,mutsig)][,.(avg=mean(n)),by=.(cluster,mutsig)][, cluster := factor(paste0('Cluster ', cluster), levels=cluster.order)],
    aes(x=mutsig, weight=avg, fill=mutsig)) +
    geom_sbs96() +
    facet_grid(cluster~., scale='free', axes='all', switch='y') +
    base.theme + manuscript +
    xlab(element_blank()) + ylab(element_blank()) +
    theme(panel.spacing=unit(1/2, unit='lines'),
        aspect.ratio=1/4,
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(), axis.text.y=element_blank(),
        strip.placement='outside',
        strip.text=element_text(margin=margin(0,0,0,0)))
#do.plot(hmsigplot, 'figure2_experimental_heatmap_cluster_sigs.pdf', width=55, height=5/4*55+5, units='mm', plot.if.exists=T)

cat('a\n')
#bottom.row <- (pg|ph|pf2i) + plot_layout(widths=c(1,4,1))
#bottom.row <- ((ph_col|pf2_heatmap) + plot_layout(widths=c(1,3))) /
#bottom.row <- (((plot_spacer()|ph_row/pf2_heatmap) + plot_layout(height=c(1, 10))) /
    #((plot_spacer()|pf2_hmsigplot) + plot_layout(widths=c(4, 1)))) # + plot_layout(height=c(2, 1))
#bottom.row <- ph_row/(pf2_heatmap|(plot_spacer()/pf2_hmsigplot) + plot_layout(width=c(4,1))) + plot_layout(height=c(1,6), guides='collect')
#bottom.row <- ph_row/((pf2_heatmap|(plot_spacer()/pf2_hmsigplot)) + plot_layout(width=c(3,1))) + plot_layout(height=c(1,8)) #, guides='collect') #& theme(legend.position='bottom')
bottom.row <- ph_row/((pf2_heatmap+pf2_hmsigplot) + plot_layout(width=c(2.25,1))) + plot_layout(height=c(1,6)) #, guides='collect') #& theme(legend.position='bottom')
cat('b\n')

fig2 <- (((pb|pc|(pd/pe)) + plot_layout(widths=c(1, 3.5, 3))) / bottom.row) +
    #plot_layout(heights=c(3,2)) &    # pre-heatmap version
    plot_layout(heights=c(2.0, 5)) &     # with heatmap. doubled bottom row height for manual move of hmsigplot
    theme(plot.margin=margin(1/2,1/2,0,1/2, unit='line'),
        #plot.background=element_rect(color='blue'),
        #panel.background=element_rect(color='red'),
        #strip.background=element_rect(color='purple'),
        #legend.position='bottom',
        plot.tag.position=c(0,1),
        text=element_text(debug=FALSE)) &
    plot_annotation(tag_levels='A', theme=theme(plot.margin=margin(0,0,0,0)))
#print(fig2)
#do.plot(fig2, 'figure2.pdf', width=180, height=105, units='mm', plot.if.exists=T)
# hack right now for heatmap is to add another row at the bottom and drag the cluster
# signatures up manually in inkscape
# so the height needs to be much larger than the final figure
do.plot(fig2, 'figure2.pdf', width=180, height=220, units='mm', plot.if.exists=T)
