varmap <- c(fdr='False discovery rate', calls='Mutation calls',  sens='Sensitivity')
# reorder by #calls in each muttype and capitalize labels
zz <- pcs[muttype != 'dnv' & variable %in% c('fdr', 'calls', 'sens'),
    .(Dataset, Tissue, Muttype=fct_shift(fct_recode(muttype, SNV='snv', Indel='indel')),
      sample=fct_reorder2(paste(Dataset, alias2), ifelse(muttype=='snv' & variable=='calls', 1, 0), value, .desc=FALSE),
      #sample=paste(Dataset, alias2),
      variable=factor(varmap[as.character(variable)], levels=varmap[c(2,3,1)]), value)]

psf2 <- ggplot(zz, aes(x=sample, y=value)) +
    facet_nested(Muttype + variable ~ Tissue, scale='free', space='free_x', switch='y', nest_line=element_line(color='black')) +
    geom_col() +
    scale_y_continuous(expand=expansion(c(0,0.05))) +
    xlab('Single cell ID') + ylab(element_blank()) +
    base.theme + manuscriptfont +
    theme(
        plot.margin=margin(1/4, 1/4, 1/4, 0, unit='line'),
        axis.text.x=element_text(angle=90, vjust=1/2, hjust=1),
        axis.title=element_text(margin=margin(0,0,0,0)),
        panel.spacing=unit(1/2, unit='line'),
        strip.placement='outside',
        strip.text=element_text(size=8, margin=margin(0.1,0.1,0.1,0.1, unit='line'))
    )
print(psf2)
do.plot(psf2, 'supp_figure2.pdf', width=7.1, height=8, plot.if.exists=TRUE)
