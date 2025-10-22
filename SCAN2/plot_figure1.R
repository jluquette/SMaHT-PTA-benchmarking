# Figure 1
# QC metrics, passing QC counts
# Include all cells regardless of QC status in this figure
pcs <- meta[per.cell.stats,,on=.(sample)]

# depth distribution that excludes bulk=0
dp.distn2 <- function(smry) {
    dptab <- decompress.dt(smry@depth.profile$dptab)

    # IMPORTANT: dptab[,-1]: the first column is where bulk depth=0. In this (and most) project(s),
    # bulk depth is 0 at generally unassayable (even unassembled, e.g. chr22's p-arm) regions. It
    # is not useful to penalize PTA for not assaying these regions, as most would interpret a
    # breadth-of-coverage number as the coverage of the assayable genome.
    ret <- data.table(sample=name(smry), bases=rowSums(dptab[,-1]), dp=as.integer(rownames(dptab)))
    ret[order(dp, decreasing=TRUE)][, cum.bases := cumsum(bases)][, cum.frac := cum.bases/max(cum.bases)][order(dp)][]
}


# MAPD example: nygc Lu_1
#pf1b <- ggplot(melt.and.annotate(mapd(sm)), aes(x=binsize, y=value, group=sample)) +
    #geom_line(linewidth=0.15, alpha=1/5) +
    #scale_x_log10(expand=c(0,0), labels=label_log()) +
    #gghighlight(n.fails<=1, use_direct_label=F, unhighlighted_params=list(colour='#888888')) +
    #ylab("Amplification non-uniformity\n(MAPD)") +
    #xlab("Bin size for MAPD calculation")
mapdz <- melt.and.annotate(mapd(sm))
p1eg_mapd <- ggplot(mapdz[n.fails < 2], aes(x=binsize, y=value)) +
    stat_summary(geom='ribbon', fun=mean,
        fun.min=function(...) quantile(..., prob=0.1),
        fun.max=function(...) quantile(..., prob=0.9), fill='grey') +
    stat_summary(geom='quantile', linewidth=1/2, col='black') +
    scale_x_log10(expand=c(0,0), labels=label_log()) +
    #coord_cartesian(xlim=c(0,1), expand=F) +
    geom_line(data=mapdz[center=='nygc' & alias=='Lu_1'], aes(x=binsize, y=value), col='firebrick2', linewidth=1/2) +
    base.theme + manuscript + #theme(aspect.ratio=1) +
    ylab('MAPD') + xlab('Bin size')


# GC example bch-broad Lu_OGQTI
#pf1c <- ggplot(melt.and.annotate(gc.bias(sm)), aes(x=gc, y=value, group=sample)) +
    #geom_line(linewidth=0.15, alpha=1/5) +
    #gghighlight(n.fails<=1, use_direct_label=F, unhighlighted_params=list(colour='#888888')) +
    #ylab("Relative read depth") + xlab("GC content (fraction bases C or G)")
gcz <- melt.and.annotate(gc.bias(sm))
p1eg_gc <- ggplot(gcz[n.fails < 2], aes(x=gc, y=value)) +
    stat_summary(geom='ribbon', fun=mean,
        fun.min=function(...) quantile(..., prob=0.1),
        fun.max=function(...) quantile(..., prob=0.9), fill='grey') +
    stat_summary(geom='quantile', linewidth=1/2, col='black') +
    #coord_cartesian(xlim=c(0,1), expand=F) +
    scale_x_continuous(expand=expansion(c(0,0))) +
    #geom_line(data=gcz[center=='bch-broad' & alias=='Lu_OGQTI'], aes(x=gc, y=value), col='firebrick2', linewidth=1/2) +
    geom_line(data=gcz[center=='bch-broad' & alias=='Lu_2'], aes(x=gc, y=value), col='firebrick2', linewidth=1/2) +
    base.theme + manuscript + #theme(aspect.ratio=1) +
    ylab('Rel. depth') + xlab('GC content')


# AB example: mayo-washu Lu_n16
#pf1d <- ggplot(melt.and.annotate(ab.distn(sm)), aes(x=af, y=value, group=sample)) +
    #geom_line(linewidth=0.15, alpha=1/5) +
    #gghighlight(n.fails<=1, use_direct_label=F, unhighlighted_params=list(colour='#888888')) +
    #ylab("Density") + xlab("Allele balance\n(VAF of het. germline SNPs)")
abz <- melt.and.annotate(ab.distn(sm))
p1eg_ab <- ggplot(abz[n.fails < 2], aes(x=af, y=value)) +
    stat_summary(geom='ribbon', fun=mean,
        fun.min=function(...) quantile(..., prob=0.1),
        fun.max=function(...) quantile(..., prob=0.9), fill='grey') +
    stat_summary(geom='quantile', linewidth=1/2, col='black') +
    coord_cartesian(xlim=c(0,1), ylim=c(0,4), expand=F) +
    geom_line(data=abz[center=='mayo-washu' & alias=='Lu_n16'], aes(x=af, y=value), col='firebrick2', linewidth=1/2) +
    base.theme + manuscript + #theme(aspect.ratio=1) +
    ylab('Density') + xlab('Allele balance')


# DP example: nygc Lu_4
dpz <- meta[rbindlist(lapply(sm, dp.distn2)),,on=.(sample)]
p1eg_dp <- ggplot(dpz[dp > 0 & n.fails < 2], aes(x=dp, y=bases/1e6)) +
    stat_summary(geom='ribbon', fun=mean,
        fun.min=function(...) quantile(..., prob=0.1),
        fun.max=function(...) quantile(..., prob=0.9), fill='grey') +
    stat_summary(geom='quantile', linewidth=1/2, col='black') +
    coord_cartesian(xlim=c(1,110), expand=F) +
    geom_line(data=dpz[dp>0 & center=='nygc' & alias=='Lu_4'], aes(x=dp, y=bases/1e6), col='firebrick2', linewidth=1/2) +
    base.theme + manuscript + #theme(aspect.ratio=1) +
    ylab('Megabases') + xlab('Seq. Depth')

pf1eg_grid <- p1eg_dp + p1eg_mapd + p1eg_ab + p1eg_gc + plot_layout(ncol=2) &
    theme(plot.margin=margin(1.5/10, 1.5/10, 1.5/10, 1.5/10, unit='line'),
        #plot.tag.position=c(0,1),
        #plot.tag.location='panel',
        #text=element_text(debug=TRUE),
        #title=element_text(debug=TRUE),
        axis.title=element_text(size=6, margin=margin(0,0,0,0))) & #, debug=TRUE)) &
    # By only providing 1 tag, we trick patchwork into not tagging all 4 panels
    plot_annotation(tag_levels=list('B'), theme=theme(plot.margin=margin(0,0,0,0)))
#do.plot(pf1eg_grid, 'figure1_b_grid.pdf', width=60, height=55, units='mm') #, plot.if.exists=F)



qcz <- meta[amp == 'PTA' & project=='smaht']
pf1qc <- ggplot(qcz, aes(x=n.fails, fill=ifelse(n.fails <= 1, 'QC pass', 'QC fail'))) +
    geom_bar() +
    base.theme + manuscriptfont +
    scale_y_continuous(expand=expansion(c(0, 0.05))) +
    geom_vline(xintercept=1.5, linewidth=0.15) +
    annotate('segment', x=0-0.4, xend=1+0.4, y=max(qcz[n.fails<=1,table(n.fails)]) + 2.5, linewidth=0.4) +
    annotate('text', x=1/2, y=max(qcz[n.fails<=1,table(n.fails)])+6, size=7, size.unit='pt',
        label=sum(qcz$n.fails<=1)) +
    annotate('segment', x=2-0.4, xend=5+0.4, y=max(qcz[n.fails>1,table(n.fails)]) + 2.5, linewidth=0.4) +
    annotate('text', x=3.5, y=max(qcz[n.fails>1,table(n.fails)]) + 6, size=7, size.unit='pt',
        label=sum(qcz$n.fails>1)) +
    guides(fill=guide_legend(position='inside', title=element_blank())) +
    theme(legend.position.inside=c(0.75,0.90), legend.key.size=unit(3.5, units='mm'), legend.key.spacing=unit(1/2, units='line')) +
    xlab("Number of failed metrics") + ylab("Single cells") +
    scale_fill_manual(values=c('firebrick2', 'black'))

# Include or exclude QC fails at this point?
#pf1depth <- ggplot(pcs[n.fails < 2 & variable=='sequencing.depth'],
pf1depth <- ggplot(pcs[variable=='sequencing.depth'],
    aes(x=fct_recode(tissue, Colon='colon', Lung='lung'), y=value)) +
    geom_boxplot(outliers=F, linewidth=0.25) +
    geom_jitter(width=0.2, size=1/4) +
    base.theme + manuscriptfont +
    xlab(element_blank()) + ylab("Mean sequencing depth per cell") #+

# not excluding failures
#pf1breadth <- ggplot(meta[rbindlist(lapply(sm, dp.distn2)),,on=.(sample)][n.fails < 2 & dp %in% c(1,10,20)],
pf1breadth <- ggplot(meta[rbindlist(lapply(sm, dp.distn2)),,on=.(sample)][dp %in% c(1,10,20)],
    aes(x=as.factor(dp), y=100*cum.frac)) +
    #aes(x=paste0('\u2265', dp), y=100*cum.frac)) +   # doesn't work because of unicode rendering PDF issue
    geom_boxplot(outliers=FALSE, linewidth=0.25) +
    geom_jitter(width=0.3, size=1/4) +
    geom_text(data=data.table(tissue=c('lung', 'colon'), cap.tissue=c('Lung', 'Colon'), x=1, y=35), aes(x=x, y=y, label=cap.tissue), size=8, size.unit='pt') +
    facet_grid(tissue~., axes='all') +
    base.theme + manuscriptfont +
    theme(panel.spacing=unit(2, units='mm'), strip.background=element_blank(), strip.text=element_blank()) +
    ylab('Percent of genome') + xlab('Sequencing depth')

# can't figure out how to size the subpatchwork grid reasonably
fig1 <- (pf1depth | pf1breadth | pf1eg_grid | pf1qc) +
    plot_layout(widths=c(2,2.5,5,3.5)) &
    theme(plot.margin=margin(1/10, 3/10, 1/10, 3/10, unit='line')) &
    # IMPORTANT!! plot_annotation(theme=) controls the margin around the entire plot.
    # can't be set any other way.
    # https://github.com/thomasp85/patchwork/issues/94
    plot_annotation(tag_levels=list(c('B', 'C', 'D', '', '', '', 'E')), theme=theme(plot.margin=margin(0,0,0,0)))
#print(fig1)
#do.plot(fig1, 'figure1.pdf', width=120, height=55, units='mm') #, plot.if.exists=F)
do.plot(fig1, 'figure1.pdf', width=180, height=55, units='mm', plot.if.exists=T)
