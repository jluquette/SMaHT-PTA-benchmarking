# must not use p[], because it has DNVs combined into single events
muts <- compute.nearest(passing(sm[meta[center=='g1' & alias=='Lu-1',sample]], combine.mnv=FALSE))
muts[, chr := factor(chr, levels=paste0('chr', c(1:22, 'X', 'Y')))]  # order chroms for the rainfall plot

pA <- plot_fancy_sbs96(muts[muttype=='snv'], mutsig) +
    base.theme + manuscriptfont +
    labs(title='Yonsei Lu_1') +
    ylab('Number of mutations') +
    theme(plot.title=element_text(hjust=1/2, size=9),
        panel.spacing=unit(0.1, units = "mm"),
        aspect.ratio=6/4,
        axis.text.x=element_text(angle = 90, size = 5, hjust = 1, vjust = 1/2, family = "mono"),
        strip.text=element_text(size=7)
    )
        #axis.line=element_line(linewidth=0.25),
        #axis.ticks=element_line(linewidth=0.25))

#pB <- ggplot(muts[chr == 'chr1' | chr=='chr2'], aes(x=pos, y=nearest, col=sbs6(mutsig))) +
pB <- ggplot(muts, aes(x=pos, y=nearest, col=sbs6(mutsig))) +
    geom_point(size=1/4) +
    facet_grid(~chr, scale='free_x', space='free_x', switch='x') +
    scale_y_log10(labels=label_log(), limits=c(1,NA)) +
    scale_x_continuous(expand=expansion(0.15)) +
    base.theme + manuscriptfont +
    xlab('Genome position') + ylab('Intermutational distance (bp)') + 
    scale_color_manual(values=sbs_cols()) +
    #guides(color=guide_legend(ncol = 2)) +
    theme(
        panel.spacing=unit(0.25, units='mm'),
        strip.placement='outside',
        strip.text=element_text(size=5, angle=90, vjust=1/2, hjust=1),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        #legend.position='inside',
        #legend.direction='horizontal',
        #legend.justification=c(0.01/2.5,0.01),
        legend.title=element_blank(),
        legend.spacing=unit(0, unit='lines'),
        legend.key.spacing=unit(1/4, unit='lines'),
        legend.key.size=unit(4, unit='mm'),           # size of the area around the point, not the point itself
        legend.box.background=element_rect(color='black'),
        legend.text=element_text(margin=margin(0,0,0,0,unit='mm')),
        legend.margin=margin(1,2,1,1,unit='mm'))  # inner margin
#print(pB)

supptab <- fread("supp_table1.csv")

pC1 <- ggplot(meta[supptab[!is.na(`Component 5`)],,on=.(sample)], aes(x=fct_reorder(.f=paste(Dataset, alias), .x=`Component 5`), weight=`Component 5`)) + geom_bar() + base.theme + manuscript + theme(axis.text.x=element_text(angle=90, vjust=1/2, hjust=1)) + xlab(element_blank()) + ylab('Component 5')

pC2 <- ggplot(meta[supptab[!is.na(`Component 5`)][p,,on=.(sample)][!is.na(`Component 5`)],,on=.(sample)], aes(x=fct_reorder(paste(Dataset, alias), `Component 5`), y=nearest, col=sbs6(mutsig))) + geom_jitter(size=1/2, stroke=0, width=0.15) + base.theme + manuscript + theme(axis.text.x=element_text(angle=90, vjust=1/2, hjust=1)) + scale_color_manual(values=sbs_cols()) + geom_hline(yintercept=32, linetype='dashed', linewidth=0.35) + xlab('Single cell ID') + ylab('Intermutational distance (bp)') + guides(color='none') + scale_y_log10(labels=label_log(), limits=c(1,NA))

pC <- pC1/pC2 + plot_layout(axes='collect_x', heights=c(1,3)) & theme(plot.margin=margin(0,0,0,0))

ps4 <- pA/pB/free(pC) + plot_layout(heights=c(1,2,4)) + plot_annotation(tag_level=list(c('A', 'B', 'C', '')), theme=theme(plot.margin=margin(0,0,0,0))) #, title='Yonsei Lu_1')
#print(ps4)
do.plot(ps4, 'supp_figure4.pdf', width=7, height=8, plot.if.exists=TRUE)
