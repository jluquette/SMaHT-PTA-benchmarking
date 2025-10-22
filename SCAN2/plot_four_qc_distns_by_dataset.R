# MAPD
p1 <- ggplot(melt.and.annotate(mapd(sm)),
             aes(x=binsize, y=value, col=index, label=alias, group=index)) +
    geom_line(data=melt.and.annotate(mapd(ns), nmeta), col=thin.col, linewidth=thin.line) +
    geom_textline(linewidth=thick.line, hjust=0.25) +
    scale_x_log10(expand=c(0,0), labels=label_log()) +
    base.theme + theme(aspect.ratio=1) +
    ylab('MAPD (amplification uniformity)') + facet_grid(center ~ tissue, scale='free') +
    guides(col='none')
do.plot(p1, 'mapd.pdf', width=12, height=24)
fwrite(melt.and.annotate(mapd(sm)), file=paste0(plot.dir, '/mapd.csv'))
fwrite(melt.and.annotate(mapd(ns), nmeta), file=paste0(plot.dir, '/mapd_neuron.csv'))


# GC bias
p2 <- ggplot(melt.and.annotate(gc.bias(sm)),
             aes(x=gc, y=value, col=index, label=alias, group=index)) +
    geom_line(data=melt.and.annotate(gc.bias(ns), nmeta), linewidth=thin.line, col=thin.col) +
    geom_textline(linewidth=thick.line) +
    base.theme + theme(aspect.ratio=1) +
    ylab('Relative amplification level') + xlab("GC content") + facet_grid(center ~ tissue, scale='free') +
    guides(col='none')
do.plot(p2, 'gc.pdf', width=12, height=24)
fwrite(melt.and.annotate(gc.bias(sm)), file=paste0(plot.dir, '/gc.csv'))
fwrite(melt.and.annotate(gc.bias(ns), nmeta), file=paste0(plot.dir, '/gc_neuron.csv'))


# Allele balance
p3 <- ggplot(melt.and.annotate(ab.distn(sm)),
             aes(x=af, y=value, col=index, label=alias, group=index)) +
    geom_line(data=melt.and.annotate(ab.distn(ns), nmeta), linewidth=thin.line, col=thin.col) +
    geom_textline(linewidth=thick.line) +
    base.theme + theme(aspect.ratio=1) +
    coord_cartesian(ylim=c(0,4)) +
    ylab('Density') + xlab("Allele balance") + facet_grid(center ~ tissue, scale='free') +
    guides(col='none')
do.plot(p3, 'ab.pdf', width=12, height=24)
fwrite(melt.and.annotate(ab.distn(sm)), file=paste0(plot.dir, '/ab.csv'))
fwrite(melt.and.annotate(ab.distn(ns), nmeta), file=paste0(plot.dir, '/ab_neuron.csv'))


# Depth and mean depth
suppressWarnings(z <- melt.and.annotate(dp.distn(sm)))
zz <- z[z[, .(dp=round(sum(dp*value)/sum(value))), by=sample],,on=.(sample,dp)]
p4 <- ggplot(z[dp>0], aes(x=dp, y=value, col=index, label=alias, group=index)) +
    geom_line(data=melt.and.annotate(dp.distn(ns), nmeta), linewidth=thin.line, col=thin.col) +
    geom_textline(hjust=1, linewidth=3/4) + geom_point(data=zz) +
    geom_text_repel(data=zz) +
    xlim(0,120) + ylim(0,1e8) +
    base.theme + theme(aspect.ratio=1) +
    ylab('Density') + xlab("Sequencing depth") + facet_grid(center ~ tissue, scale='free') +
    guides(col='none')
do.plot(p4, 'dp.pdf', width=12, height=24)
fwrite(melt.and.annotate(dp.distn(sm)), file=paste0(plot.dir, '/dp.csv'))
fwrite(melt.and.annotate(dp.distn(ns), nmeta), file=paste0(plot.dir, '/dp_neuron.csv'))


# "reason" matches the codes in the QC fail file
cz <- rbind(
    melt.and.annotate(mapd(sm))[, c('metric', 'reason') := list('MAPD', 'MAPD')][, .(sample, Dataset, alias, reason, metric, x=binsize, y=value)],
    melt.and.annotate(gc.bias(sm))[, c('metric', 'reason') := list('GC bias', 'GC')][, .(sample, Dataset, alias, reason, metric, x=gc, y=value)],
    melt.and.annotate(ab.distn(sm))[, c('metric', 'reason') := list('Allele balance', 'AB')][, .(sample, Dataset, alias, reason, metric, x=af, y=value)],
    melt.and.annotate(dp.distn(sm))[, c('metric', 'reason') := list('Sequencing depth', 'DP')][, .(sample, Dataset, alias, reason, metric, x=dp, y=value)][x > 0 & x <= 120])

# annotate with failures
cz <- qcreasons[, fail := TRUE][cz,,on=.(sample,reason)][is.na(fail), fail := FALSE]

pc <- ggplot(cz, aes(x=x, y=y, col=fail, label=alias, group=sample)) +
    geom_textline(linewidth=3/4) +
    base.theme + manuscriptfont +
    theme(aspect.ratio=1, panel.spacing=unit(1/2, unit='lines'), legend.position='bottom', plot.margins=margin(0,0,0,0)) +
    facet_grid2(Dataset ~ metric, scale='free', independent='all') +
    xlab(element_blank())+ ylab(element_blank())
do.plot(pc, 'supp_figure1.pdf', width=7.1, height=8.5, plot.if.exists=TRUE)
