# "reason" matches the codes in the QC fail file
cz <- rbind(
    melt.and.annotate(mapd(sm))[, c('metric', 'reason') := list('MAPD', 'MAPD')][, .(sample, Dataset, alias, reason, metric, x=log10(binsize), y=value)],
    melt.and.annotate(gc.bias(sm))[, c('metric', 'reason') := list('GC bias', 'GC')][, .(sample, Dataset, alias, reason, metric, x=gc, y=value)],
    melt.and.annotate(ab.distn(sm))[, c('metric', 'reason') := list('Allele balance', 'AB')][, .(sample, Dataset, alias, reason, metric, x=af, y=pmin(5,value))],
    melt.and.annotate(dp.distn(sm))[, c('metric', 'reason') := list('Sequencing depth', 'DP')][, .(sample, Dataset, alias, reason, metric, x=dp, y=value)][x > 0 & x <= 120])

cz[Dataset=="Yonsei", alias := gsub(pattern='-', replace='_', alias)]

# annotate with failures
cz <- qcreasons[, fail := TRUE][cz,,on=.(up_id_header=sample,alias,reason)][is.na(fail), fail := FALSE][, sample := up_id_header][]

pc <- ggplot(cz, aes(x=x, y=y, label=alias, group=sample)) +
    geom_textline(linewidth=1/3, size=2) +
    gghighlight(paste(sample, reason) %in% cz[fail==TRUE, paste(sample, reason)], calculate_per_facet=TRUE) +
    base.theme + manuscriptfont +
    theme(aspect.ratio=1,
        strip.placement='outside',
        panel.spacing=unit(1/2, unit='lines'),
        #legend.position='bottom',
        plot.margins=margin(0,0,0,0)) +
    facet_grid2(metric ~ Dataset, scale='free', independent='x', switch='y') +
    xlab(element_blank())+ ylab(element_blank())
do.plot(pc, 'supp_figure1.pdf', width=7.0, height=6.2, plot.if.exists=TRUE)
