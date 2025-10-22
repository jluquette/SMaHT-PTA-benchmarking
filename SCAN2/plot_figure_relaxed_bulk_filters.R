# Compare relaxed bulk filtered calls to standard calling
# Heuristic: assume all relaxed calls are FPs to put a ceiling on FDR
z <- meta[rbind(p,pr),,on=.(sample)]
zz <- meta[dcast(z[, .N, by=.(analysis, sample, muttype)], sample+muttype ~ analysis, value.var='N')[, fdr := 1-standard/relaxed],,on=.(sample)]

p1 <- ggplot(zz, aes(standard, relaxed, col=tissue)) +
    geom_point(size=1/3) +
    facet_wrap(~muttype, scale='free') +
    geom_abline()

p2 <- ggplot(zz, aes(x=sample, y=fdr, fill=tissue)) +
    geom_col() +
    facet_grid(muttype~., scale='free') +
    ylab('Worst-case false discovery rate') +
    theme(axis.text.x=element_text(angle=90, vjust=1/2, hjust=1, size=5))

p3 <- ggplot(zz, aes(x=interaction(tissue,muttype, sep='\n'), y=fdr, col=tissue)) +
    geom_boxplot(outliers=FALSE) +
    geom_jitter(size=1/3, width=0.2) +
    ylab('Worst-case false discovery rate') + xlab(element_blank())

rbplot <- (p1/((p2|p3) + plot_layout(width=c(4,1)))) +
    plot_annotation(title='Standard SCAN2 vs. relaxed bulk filters to allow clonal variant calling',
        subtitle='Worst-case FDR: assume ALL new calls from relaxed filters are false positives') +
    plot_layout(guides='collect') & manuscriptfont
do.plot(rbplot, 'effect_of_relaxed_bulk_filters.pdf', width=180, height=140, units='mm', plot.if.exists=F)



# Might be a good place for violin plots. The <100bp clusters are not as visible now that
# width=0.1 to make the boxplot edges more visible.
# must reorder here because only want to consider standard calling
z2 <- copy(z)[, sample := fct_reorder(sample, ifelse(analysis=='standard',nearest,Inf), .fun=function(x) sum(x<500 & x>1))]
imd.plot <- ggplot(z2[muttype=='snv'], aes(x=sample, y=nearest)) +
    geom_boxplot(outliers=FALSE, linewidth=0.2) +
    geom_jitter(aes(col=substr(mutsig, 5, 7)), width=0.20, stroke=0, size=1/5) +
    facet_grid(analysis~tissue, scale='free_x', space='free_x') +
    scale_y_log10(labels=label_log()) +
    scale_color_manual(values=sbs_cols()) +
    guides(color=guide_legend(title='Base change', override.aes=list(size=2), nrow=1)) +
    xlab('Sample') + ylab('Intermutational distance (IMD)') +
    manuscript + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2, size=5), legend.position='bottom')
do.plot(imd.plot, 'intermutation_distance.relaxed_vs_standard.snv_only.pdf', width=180, height=120, units='mm', plot.if.exists=F)
