z <- fread('hdp_from_tim/release_2/hdp_sigs_combined_v2.txt')
zz <- rbindlist(lapply(3:8, function(i)
    data.table(component=paste('Component', i-2),
               muttype=factor(rep(c('SNV', 'Indel', 'DNV'), times=c(96, 83, 78)),
               levels=c('SNV', 'Indel', 'DNV')),
               mutsig=cmp257(levels(cmp257())), value=z[(1:258)[-(96+83+1)]][[i]])
))

ps3 <- ggplot(zz, aes(x=mutsig, fill=mutsig, weight=value*100)) +
    geom_bar() + 
    scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
    scale_y_continuous(expand=expansion(c(0,0.05))) +
    facet_grid2(component ~ muttype, scale='free', switch='both', space='free_x', independent='y') +
    #ylab("Percent of mutations") + xlab(element_blank()) +
    ylab(element_blank()) + xlab(element_blank()) +
    base.theme + manuscriptfont +
    guides(x=guide_axis(n.dodge=2)) +
    theme( #aspect.ratio=1/2.0,
        plot.margin=margin(0.25,0.25,0.25,0.05, unit='line'),
        strip.text.y=element_text(size=8, face='bold', margin=margin(0,0,0,0)),
        strip.text.x=element_text(size=8, face='bold'),
        strip.placement='outside',
        panel.spacing.x=unit(0.5, unit='line'),
        panel.spacing.y=unit(0.5, unit='line'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_text(size=4.5, family='mono', angle=90, vjust=1/2, hjust=1),
        axis.text.y=element_text(size=6, angle=90, vjust=1, hjust=1/2, margin=margin(0,0,0,0)) )
do.plot(ps3, 'supp_figure3.pdf', width=7.1, height=6.5, plot.if.exists=TRUE)
