library(ggplot2)
library(patchwork)
library(geomtextpath)
library(scales)
library(forcats)
library(ggrepel)
library(scan2)
library(simplesigs)
library(ggforce)
library(concaveman)
library(mclust)
library(umap)
library(gghighlight)
library(ggbreak)
library(ggcorrheatmap)
library(ggh4x)

tissue.colors <- c(colon='#9C2525', Colon='#9C2525',
    lung='#7BC8EA', Lung='#7BC8EA')

centers=c('bch-broad', 'yale-bcm', 'g1', 'mayo-washu', 'nygc')
center.to.dataset.map <- c(
    `bch-broad`='BCH-Broad',
    `yale-bcm`='Yale-BCM',
    `g1`='Yonsei',
    `mayo-washu`='Mayo-WashU',
    `nygc`='NYGC')

center.colors <- c(
    `bch-broad`='#FCD201',
    `BCH-Broad`='#FCD201',
    `yale-bcm`='#20823B',
    `Yale-BCM`='#20823B',
    `g1`='#BC8F61',
    `Yonsei`='#BC8F61',
    `mayo-washu`='#427CBF',
    `Mayo-WashU`='#427CBF',
    `nygc`='#764897',
    `NYGC`='#764897'
)

plot.dir <- 'plots/analysis_september_10_2025'
replot <- FALSE # If a plot already exists, replot it?
#replot <- TRUE   # If a plot already exists, replot it?

# Set to true to load all summary objects and metadata EVERY TIME
# the script is run. This is not practical for developing the plots.
reload.data <- FALSE
#reload.data <- TRUE


# Sizes suitable for manuscript figures
thin.line <- 0.25
thick.line <- 0.75
thin.col <- 'lightgrey'

manuscriptfont <- theme(
    axis.title=element_text(size=7),
    axis.text=element_text(size=6),
    strip.text=element_text(size=8),
    legend.title=element_text(size=8),
    legend.text=element_text(size=7),
    plot.tag=element_text(size=9, face='bold')
)
manuscriptlines <- theme(
    line=element_line(linewidth=0.15)
)
manuscript <- manuscriptfont + manuscriptlines


# Sizes suitable for powerpoint slides
base.text=20
large.text=22
signif.size=1.6 # size of the p-value labels on significance tests. not sure what the unit is

smallerfont <- theme(
    axis.title=element_text(size=large.text-8),
    axis.text=element_text(size=base.text-8),
    strip.text=element_text(size=large.text-8),
    legend.title=element_text(size=large.text-8),
    legend.text=element_text(size=base.text-8)
)

bigfont <- theme(
    axis.title=element_text(size=large.text),
    axis.text=element_text(size=base.text),
    strip.text=element_text(size=large.text),
    legend.title=element_text(size=large.text),
    legend.text=element_text(size=base.text)
)

base.theme <- list(
    theme_classic(),
    bigfont,
    theme(axis.line=element_line(linewidth=0.25),
        axis.ticks=element_line(linewidth=0.25),
        strip.background=element_blank(),
        panel.spacing=unit(2.25, 'lines'))
)


if (TRUE | !'meta' %in% ls()) {
    meta <- fread("sample_metadata.csv")
    meta[, Tissue := ifelse(tissue=='lung', 'Lung', ifelse(tissue=='colon', 'Colon', tissue))]
    meta[, Dataset := center.to.dataset.map[center]]
    #meta[, bam := NULL]
    # index: a much better option for coloring than alias because the same
    #   index value is shared across all centers/tissues while aliases may
    #   be slightly different. This allows color sharing across facets.
    # If converted to a factor, easier to have ggplot color/sort the output
    #meta[, index := as.factor(1:.N), by=.(tissue,center)]
    # short alias for manuscript figures
    meta[, alias2 := paste0(ifelse(tissue=='lung', 'Lu', 'Co'), '_', 1:.N), by=tissue]


    qcreasons <- fread('qc_fails_with_reasons.txt')[, sample := up_id_header]
    # Fails more than one visual inspection
    qc.fail <- qcreasons[,.N,by=.(sample)][N>1,sample]
    meta <- merge(meta, qcreasons[, .(n.fails=.N), by=.(sample)], by='sample', all.x=T)[is.na(n.fails), n.fails:=0]
    
    nmeta <- fread("neuron_metadata.csv") # for Luquette et al 2022 neuron comparison
    nmeta[, alias := paste0('n', 1:.N)]
    nmeta[, index := as.factor(1:.N)]
}


if (reload.data | !all(c('sm.std', 'ns') %in% ls())) {
#if (reload.data | !all(c('sm.std') %in% ls())) {
    # Standard calling (i.e., no relaxing bulk filters)
    #std.summary.dir <- 'objs/summary_mock_rescue'
    std.summary.dir <- 'objs/summary'
    std.summary.files <- list.files(path=std.summary.dir, pattern='.*.rda', full.names=TRUE)
    std.summary.names <- gsub('.rda', '', list.files(path=std.summary.dir, pattern='.*.rda'))

    # Relaxed bulk filters to improve clonal/truth set variant recovery
    #rb.summary.dir <- 'objs/summary_mock_rescue_bulk_binom_prob_1e-5'
    rb.summary.dir <- 'objs/summary_bulk_binom_prob_1e-5'
    rb.summary.files <- list.files(path=rb.summary.dir, pattern='.*.rda', full.names=TRUE)

    # 52 PTA neurons from Luquette et al 2022. For comparison of amplification quality
    n.summary.dir <- '/n/data1/hms/dbmi/park/jluquette/SMaHT/pta_pilot/objs/neurons_PTA_2022'
    n.summary.files <- list.files(path=n.summary.dir, pattern='.*.rda', full.names=TRUE) #[1:2] # hack for now - just read in 2 neurons to save loading time

    sm.std <- load.summary(paths=std.summary.files)
    # Not useful right now
    #sm.rb <- load.summary(paths=rb.summary.files)
    ns <- load.summary(paths=n.summary.files)

    sm <- sm.std

    # These tables are created once in run_digest.sh, we do not rely on the
    # summary objects for calls.
    p <- fread("all_SCAN2_pass_snv_indel_dnv.csv.gz")[muttype != 'mnv']
    p[, analysis := 'standard']
    compute.nearest(p)
    pr <- fread("all_SCAN2_pass_snv_indel_dnv.relaxed.csv.gz")[muttype != 'mnv']
    pr[, analysis := 'relaxed']
    compute.nearest(pr)

    data.loaded <- TRUE
    print(gc())
}


cat("WARNING: per.cell.stats: for calculating FPRs/FDRs, using static table of mutations from scan2 rescue with rescued sites removed, pre-digest_calls.R filtering\n")
per.cell.stats <- callstats(sm, pdata=p)
#neuron.per.cell.stats <- rbindlist(lapply(ns, stats.one, p=passing(ns)))

just.melt <- function(x, ...) melt(data.table(x), id.vars=1, variable.name='sample', ...)
# this version probably works for more SCAN2 accessors, but don't want to
# break existing code.
just.melt2 <- function(x, ...) melt(data.table(x, keep.rownames='sample'), id.vars=1, ...)

# annotate with the new single cell metadata by default. for neurons use nmeta
melt.and.annotate <- function(x, metadata=meta, ...) metadata[just.melt(x), , on=.(sample)]

do.plot <- function(p, file, dev=pdf, plot.if.exists=replot, ...) {
    file <- paste0(plot.dir, '/', file)
    exists <- file.exists(file)
    if (!exists | plot.if.exists == TRUE) {
        if (exists)
            cat(paste0('replot=', plot.if.exists, ': overwriting existing plot ', file, '\n'))
        ggsave(p, dev=dev, dpi=600, file=file, ...)
    } else {
        cat(paste0('replot=', plot.if.exists, ': skipping existing plot ', file, '\n'))
    }
}


cat('start plot_figure1.R\n')
source('plot_figure1.R')
#stop('figure 1 done')

cat('start plot_figure2.R\n')
source('plot_figure2.R')
stop('figure 2 done')

cat('start plot_figure_relaxed_bulk_filters.R\n')
source('plot_figure_relaxed_bulk_filters.R')
#stop('figure relaxed bulk filters done')

cat('start plot_four_qc_distns_by_dataset.R\n')
source('plot_four_qc_distns_by_dataset.R')
#stop('figure four QC distns (by dataset) done')


# CNVs
z <- meta[melt(data.table(binned.counts(sm, type='cnv'), keep.rownames='chr')[, abs.pos := 1:.N], id.vars=c('chr', 'pos', 'abs.pos'), variable.name='sample'),,on=.(sample)][, value := pmin(value, 6)]

labels <- z[sample==sample[1],round(mean(abs.pos)),by=chr][seq(1,.N,2), setNames(chr, V1)]

for (cn in centers) {
    for (ttype in c('colon', 'lung')) {
        nrows=5
        if (cn=='yale-bcm' & ttype=='colon')  # 18 cells rather than <=10
            nrows=6
        if (cn=='mayo-washu' & ttype=='lung')  # 26 cells rather than <=10
            nrows=9
        if (nrow(z[tissue==ttype & center==cn]) > 0) {
            p5 <- ggplot(z[tissue==ttype & center==cn], aes(x=abs.pos, y=value, col=chr)) +
                geom_point(size=1/20) +
                facet_wrap(~alias, nrow=nrows, strip.position='left', axes='all', axis.labels='all_y') +
                base.theme + smallerfont +
                theme(aspect.ratio=1/7, legend.position='none', axis.text.x=element_text(angle=90, hjust=1, vjust=1/2),
                    strip.placement='outside', strip.text=element_text(size=base.text-8), panel.spacing=unit(1/4, units='lines')) +
                scale_color_manual(values=setNames(rep(c('black', 'gray50'), 12), paste0('chr', c(1:22,'X', 'Y')))) +
                scale_x_continuous(breaks=as.integer(names(labels)), labels=labels, expand=c(0,0)) +
                xlab(element_blank()) + ylab('Copy number')
            do.plot(p5, dev=png, file=paste0('cnv_', cn, '_', ttype, '.png'), width=12.5, height=6, res=300)
        } else {
            cat("skipping CNV plot for non-existant pair", cn, ttype, '\n')
        }
    }
}



# Counts/burdens of SNVs/indels, useful for several plots
counts <- p[,.(snv.calls=sum(muttype=='snv'), indel.calls=sum(muttype=='indel'), dbs.calls=sum(muttype=='dbs')),by=.(sample)]
burdens <- data.table(mutburden(sm), keep.rownames='sample')[, .(sample=sample, snv.burden=snv, indel.burden=indel)]
z <- meta[melt(counts[burdens,,on=.(sample)], id.vars=1, measure.vars=measure(muttype, type, sep='.')),,on=.(sample)]
fwrite(z, file=paste0(plot.dir, '/muts_burdens.csv'))

# Counts/burdens of somatic SNVs (NOT with relaxed bulk)
# bch-broad sample=OGQTH really low (= doublet, incompatible with mut/art mdoel?)
p6 <- ggplot(z[muttype=='snv'], aes(x=center, y=value, group=center)) +
    geom_boxplot(outliers=FALSE) +
    geom_point(aes(col=index)) + geom_text_repel(aes(col=index, label=alias)) +  # labels can't follow jitter
    facet_wrap(tissue~type, scale='free') +
    base.theme + xlab(element_blank()) + ylab("Somatic SNVs") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2), aspect.ratio=1) +
    guides(col='none')
do.plot(p6, 'somatic_SNV_counts_and_burdens.pdf', width=12, height=13, plot.if.exists=F)

# Counts/burdens of somatic indels
p7 <- ggplot(z[muttype=='indel'], aes(x=center, y=value, group=center)) +
    geom_boxplot(outliers=FALSE) + geom_point(aes(col=index)) +
    geom_text_repel(aes(col=index, label=alias)) +
    facet_wrap(tissue~type, scale='free') +
    base.theme + xlab(element_blank()) + ylab("Somatic Indels") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=1/2), aspect.ratio=1) +
    guides(col='none')
do.plot(p7, 'somatic_indel_counts_and_burdens.pdf', width=12, height=13, plot.if.exists=F)


# SNV calls vs. indel calls
p7a <- ggplot(dcast(z, ... ~ muttype + type, value.var='value', sep='.'),
    aes(x=snv.calls, y=indel.calls, col=index, label=alias)) +
    geom_point() +
    geom_text_repel() +
    facet_grid(tissue~center) +
    guides(col='none') +
    base.theme + xlab("Somatic SNVs called") + ylab("Somatic indels called") +
    theme(aspect.ratio=1)
do.plot(p7a, 'somatic_snv_vs_indel_counts.pdf', width=12, height=6)

# SNV burdens vs. indel burdens
p7b <- ggplot(dcast(z, ... ~ muttype + type, value.var='value', sep='.'),
    aes(x=snv.burden, y=indel.burden, col=index, label=alias)) +
    geom_point() +
    geom_text_repel() +
    facet_grid(tissue~center) +
    guides(col='none') +
    base.theme + xlab("Somatic SNV burden") + ylab("Somatic indel burden") +
    theme(aspect.ratio=1)
do.plot(p7b, 'somatic_snv_vs_indel_burdens.pdf', width=12, height=6)



# AF distns of called somatic mutations
# is OGQTH a doublet?
p8 <- ggplot(meta[p,,on=.(sample)], aes(x=pmax(af,1-af), col=index, label=alias, group=index)) +
    geom_textdensity(linewidth=thick.line, hjust=0.2) +
    facet_grid(center ~ tissue, scale='free') +
    base.theme +
    xlab('Mirrored VAF (>1/2)') + ylab('Density') + theme(aspect.ratio=1) +
    guides(col='none')
do.plot(p8, 'somatic_vaf.pdf', width=12, height=24) #, plot.if.exists=TRUE)



#############################################################################
# signatures plotted separately per type (SBS, ID, DBS)
#   OPTIONALLY: with transcribed strand - requires SigProfilerMatrixGenerator
#############################################################################

# weight signatures in x by sample and by center
weight <- function(x)
    melt(x[, weight.sample := 1/.N, by=.(sample)][, weight.center := length(unique(sample))/.N, by=center],
        variable.name='weighting', value.name='w',
        measure.vars=c('weight.sample', 'weight.center'))

snv <- weight(meta[p[muttype=='snv'][, mutsig := sbs96(mutsig)],,on=.(sample)])
id <- weight(meta[p[muttype=='indel'][, mutsig := id83(mutsig)],,on=.(sample)])
dbs <- weight(meta[p[muttype=='dbs'][, mutsig := dbs78(mutsig)],,on=.(sample)])

do.tx <- FALSE
if (do.tx) {
    # Read larger classifications from SigProfilerMatrixGenerator
    stx <- fread('mutsigs/sbs.txt', col.names=c('sample', 'chr', 'pos', 'ctx'))[, strand := tx(ctx)][]
    itx <- fread('mutsigs/id.txt', col.names=c('sample', 'chr', 'pos', 'ctx'))[, strand := tx(ctx)][]
    dtx <- fread('mutsigs/dbs.txt', col.names=c('sample', 'chr', 'pos', 'ctx'))[, strand := tx(ctx)][]
    # deleted the joins, need to join on=.(sample, chr, pos)
}

if (FALSE) {
for (tt in c('lung', 'colon')) {
    p9a <- ggplot(snv[tissue==tt], aes(x=mutsig, weight=w, fill=mutsig)) +
        facet_grid(center~weighting, scale='free_y') +
        base.theme + smallerfont +
        geom_sbs96()
    do.plot(p9a, paste0('signature_by_center_snv_', tt, '.pdf'), width=12, height=8, plot.if.exists=F)

    if (do.tx) {
        p9a2 <- ggplot(snv[tissue==tt][, mutsig := tx(strand, mutsig, plot=T)][!is.na(mutsig)],
            aes(x=mutsig, weight=w, fill=mutsig)) +
            facet_grid(center~weighting, scale='free_y') +
            base.theme + smallerfont +
            geom_sbs96(tx=TRUE)
        do.plot(p9a2, paste0('signature_by_center_snv_tx_', tt, '.pdf'), width=12, height=8, plot.if.exists=F)
    }


    # Indels
    p9c <- ggplot(id[tissue==tt], aes(x=mutsig, weight=w, fill=mutsig)) +
        facet_grid(center~weighting, scale='free_y') +
        base.theme + smallerfont +
        geom_id83()
    do.plot(p9c, paste0('signature_by_center_indel_', tt, '.pdf'), width=12, height=8, plot.if.exists=F)

    if (do.tx) {
        p9c2 <- ggplot(id[tissue==tt][, mutsig := tx(strand, mutsig, plot=TRUE)][!is.na(mutsig)],
            aes(x=mutsig, weight=w, fill=mutsig)) +
            facet_grid(center~weighting, scale='free_y') +
            base.theme + smallerfont +
            geom_id83(tx=TRUE)
        do.plot(p9c2, paste0('signature_by_center_indel_tx_', tt, '.pdf'), width=12, height=8, plot.if.exists=F)
    }


    # DINUCs
    # TRY: points + barplots where points show per-cell values against barplot average
    #   - this is good, but to do it the geom_bar(weight) strategy in simplesigs needs
    #     to be replaced by a new strategy. first, the spectrum of each sample needs to
    #     be computed as a fraction, then the bar height should be the mean fraction across
    #     samples, not the stacked bar of weights.
    p9e <- ggplot(dbs[tissue==tt], aes(x=mutsig, weight=w, fill=substr(mutsig, 1, 2))) +
        facet_grid(center ~ weighting, scale='free_y', axes='all') +
        base.theme + smallerfont +
        geom_dbs78(guide=TRUE) + theme(axis.text.x=element_text(size=rel(0.4)))
    do.plot(p9e, paste0('signature_by_center_dbs_', tt, '.pdf'), width=12, height=8, plot.if.exists=F)

    if (do.tx) {
        p9e2 <- ggplot(dbs[tissue==tt][, mutsig := tx(strand, mutsig, plot=TRUE)][!is.na(mutsig)],
            aes(x=mutsig, weight=w, fill=mutsig)) +
            facet_grid(center ~ weighting, scale='free_y', axes='all') +
            base.theme + smallerfont +
            geom_dbs78(guide=TRUE, tx=TRUE) + theme(axis.text.x=element_text(size=rel(0.4)))
        do.plot(p9e2, paste0('signature_by_center_dbs_tx_', tt, '.pdf'), width=12, height=8, plot.if.exists=F)
    }
}

stop('signature plots')
}

#############################################################################
# COMPOSITE signatures
#   OPTIONALLY: with transcribed strand - requires SigProfilerMatrixGenerator
#############################################################################

# reweight such that each signature type sums to 1 so that singleton DBSes don't
# blow up the y-scale
#reform <- function(x) {
    #y <- x[weighting == 'weight.sample'][, .(sample, tissue, center, strand, mutsig, muttype)]
    #y <- dcast(y, sample + tissue + center ~ muttype + mutsig, drop=TRUE, fun.aggregate=length)
    #y <- melt(y, id.vars=c('sample', 'tissue', 'center'), value.name='w', measure.vars=measure(muttype, mutsig, sep='_'))
    #y[, w := as.numeric(w)][]
#}

#allsigs <- meta[reform(rbind(snv, id, dbs))[, w := w/ifelse(w==0, 1, max(w)), by=.(sample, muttype)],,on=.(sample,center,tissue)]
#allsigs[, mutsig := factor(mutsig, levels=names(c(sbs96_cols(), id83_cols(), dbs78_cols())))]
#allsigs[, mutsig := cmp257(mutsig)][, muttype := factor(muttype, levels=c('snv', 'indel','dbs'))]

# This is critical for plotting composite signatures.
# It produces count=0 entries for every factor level in cmp257, which prevents
# dropping bars of 0 count.
mutsig_to_muttype <- c(
    setNames(rep('snv', 96), levels(sbs96())), #'snv'),
    setNames(rep('indel', 83), levels(id83())), #'indel'),
    setNames(rep('dbs', 78), levels(dbs78())) #, 'dbs')
)
allsigs <- melt(copy(p)[, mutsig := cmp257(mutsig)][, as.list(table(mutsig)), by=.(sample)], id.vars='sample', variable.name='mutsig', value.name='count')
allsigs <- meta[allsigs,, on=.(sample)][, muttype := factor(mutsig_to_muttype[mutsig], levels=c('snv', 'indel','dbs'))]


if (do.tx) {
    allsigs.tx <- meta[reform(rbind(snv, id, dbs)[, mutsig := tx(strand, mutsig, plot=TRUE)][!is.na(mutsig)])[, w := w/ifelse(w==0, 1, max(w)), by=.(sample, muttype)],,on=.(sample,center,tissue)]
    # order the signature channels
    allsigs.tx[, mutsig := factor(mutsig, levels=apply(expand.grid(c('T','U'), names(c(sbs96_cols(), id83_cols(), dbs78_cols()))),1,paste,collapse=':'))]
}


for (tt in c('lung', 'colon')) {
    for (cn in centers) {
        if (cn == 'g1' & tt == 'colon')
            next

        # geom_cmp257 doesn't quite work, maybe because of scale_x_discrete(drop=F)?
        #p10 <- ggplot(allsigs[center==cn & tissue == tt], aes(x=mutsig,weight=w, fill=mutsig)) +
        p10 <- ggplot(allsigs[center==cn & tissue == tt], aes(x=mutsig, weight=count, fill=mutsig)) +
            geom_bar() +
            facet_wrap(~interaction(muttype, alias), scale='free', ncol=3) +
            scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
            ylab(element_blank()) + xlab(element_blank()) +
            theme(aspect.ratio=1/4) + smallerfont +
            scale_y_continuous(expand=expansion(c(0,0.05))) +
            theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            labs(title=cn)
            #facet_wrap(factor(alias, levels=sample.order.by.dbs2)~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free', ncol=3, drop=F) +
            #theme(axis.text.x=element_text(size=rel(0.4), angle=90, vjust=1/2, hjust=1)) +
        do.plot(p10, paste0('composite_signatures_by_sample_', cn, '_', tt, '.pdf'), width=14,
            # scale height with #samples
            height=2 + length(unique(meta[center==cn & tissue==tt,sample])),
            plot.if.exists=F)

        if (do.tx) {
            stop('not updated to use raw counts rather than weights')
            # tx strand bias
            p10b <- ggplot(allsigs.tx[center==cn & tissue == tt], aes(x=mutsig, weight=w, fill=mutsig)) +
                geom_bar() +
                facet_grid(factor(alias, levels=sample.order.by.dbs2)~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
                scale_fill_manual(values=tx_cols(c(sbs96_cols(), id83_cols(), dbs78_cols())), guide='none') +
                ylab(element_blank()) + xlab(element_blank()) +
                theme(aspect.ratio=1/4) + smallerfont +
                scale_y_continuous(expand=expansion(c(0,0.05))) +
                theme(axis.text.x=element_text(size=rel(0.4), angle=90, vjust=1/2, hjust=1)) +
                labs(title=cn)
            do.plot(p10b, paste0('composite_signatures_by_sample_', cn, '_', tt, '_tx.pdf'), width=14,
                # scale height with #samples
                height=2 + length(unique(meta[center==cn & tissue==tt,sample])),
                plot.if.exists=F)
        }
    }

    #p10c <- ggplot(allsigs[tissue == tt][, w := w/max(w), by=.(center,muttype)],
    p10c <- ggplot(allsigs[tissue == tt],
        aes(x=mutsig, weight=count, fill=mutsig)) +
        geom_bar() +
        facet_wrap(~interaction(muttype, center), scale='free', ncol=3) +
        scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
        ylab(element_blank()) + xlab(element_blank()) +
        theme(aspect.ratio=1/4) + smallerfont +
        scale_y_continuous(expand=expansion(c(0,0.05))) +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #theme(axis.text.x=element_text(size=rel(0.4), angle=90, vjust=1/2, hjust=1))
        #facet_grid(center~factor(muttype, levels=c('snv', 'indel', 'dbs')), scale='free') +
    do.plot(p10c, paste0('composite_signatures_by_center_', tt, '.pdf'), width=9, height=7,
        plot.if.exists=F)

    if (do.tx) {
        # tx strand bias
        stop('not updated to use raw counts rather than weights')
        p10d <- ggplot(allsigs.tx[tissue == tt][, w := w/max(w), by=.(center,muttype)],
            aes(x=mutsig, weight=w, fill=mutsig)) +
            geom_bar() +
            facet_grid(center~factor(muttype, levels=c('snv', 'indel', 'dbs')), scale='free') +
            scale_fill_manual(values=tx_cols(c(sbs96_cols(), id83_cols(), dbs78_cols())), guide='none') +
            ylab(element_blank()) + xlab(element_blank()) +
            theme(aspect.ratio=1/4) + smallerfont +
            scale_y_continuous(expand=expansion(c(0,0.05))) +
            theme(axis.text.x=element_text(size=rel(0.4), angle=90, vjust=1/2, hjust=1))
        do.plot(p10d, paste0('composite_signatures_by_center_', tt, '_tx.pdf'), width=9, height=7,
            plot.if.exists=F)
    }
}



# FROM TIM - HDP signature weights
# not used anywhere currently, check again later
if (FALSE) {
ft <- melt(fread("/n/data1/hms/dbmi/park/jluquette/SMaHT/pta_pilot/from_tim_mean_assignment_hdp_combined_sc.txt"), id.vars=c(1), variable.name='sample')[, comp:=V1][,V1 := NULL]
ft[, center := gsub('.', '-', sapply(strsplit(as.character(sample), '_', fixed=T),`[`,1), fixed=T)][][, alias := sapply(strsplit(as.character(sample), '_', fixed=T), function(x) paste(x[-1], collapse='_'))]
ft2 <- dcast(ft, center+alias ~ comp, value.var='value')
}


#allsigs2 <- allsigs[!(paste(center, alias) %in% qc.fail)]

# Cluster by spectra
if (FALSE|T) {
# sample not in QC
    m <- as.matrix(dcast(allsigs[T|n.fails < 2], sample~mutsig, value.var='count'), rownames='sample')
    #m[is.na(m)]<-0
    #u <- as.data.table(umap(m)$layout, keep.rownames='sample')
    # does getting rid of mostly-0 channels do anything useful?
    u <- as.data.table(umap(m)$layout, keep.rownames='sample')
    bic <- mclustBIC(u[,.(V1,V2)])
    mod <- Mclust(u[,.(V1,V2)], x=bic)
    u[, class := mod$classification]
    #z <- meta[u,,on=.(sample)]
    # no difference w.r.t. allsigs, but metadata has been updated with cluster class
    allsigs <- u[allsigs,,on=.(sample)]
    if (do.tx)
        allsigs.tx <- u[allsigs.tx,,on=.(sample)]

    # "Genetic UMAP"
    p11c <- ggplot(allsigs, aes(x=V1, y=V2, shape=tissue)) +
        geom_mark_hull(aes(group=class, label=class), concavity=10, fill='lightgrey', alpha=0.25) +
        geom_point(aes(col=n.fails>1), size=3) +
        base.theme + smallerfont
    do.plot(p11c, 'genetic_umap_allsigs.pdf', width=9, height=9, plot.if.exists=TRUE)

    # Signatures per cluster
    p11d <- ggplot(allsigs, aes(x=mutsig, weight=count, fill=mutsig)) +
        geom_bar() +
        scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
        facet_wrap(~interaction(muttype, class), scale='free', ncol=3) +
        theme(aspect.ratio=1/4, axis.text.x=element_blank(), axis.ticks.x=element_blank())
        #facet_grid(class~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
    do.plot(p11d, 'genetic_umap_signatures_per_cluster.pdf', width=12, height=8, plot.if.exists=TRUE)

    # don't remember ever looking at these transposed plots
    if (FALSE) {
        p11d <- ggplot(allsigs, aes(x=mutsig, weight=count, fill=mutsig)) +
            geom_bar() +
            facet_grid(factor(muttype, levels=c('snv', 'indel','dbs')) ~ class, scale='free', space='free') +
            scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') #+
            #theme(aspect.ratio=1/4)
        do.plot(p11d, 'genetic_umap_signatures_per_cluster.t.pdf', width=36, height=6.5, plot.if.exists=TRUE)
    }

    if (do.tx) {
        p11d2 <- ggplot(allsigs.tx, aes(x=mutsig, weight=w, fill=mutsig)) +
            geom_bar() +
            facet_grid(class~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
            scale_fill_manual(values=tx_cols(c(sbs96_cols(), id83_cols(), dbs78_cols())), guide='none') +
            theme(aspect.ratio=1/4)
        do.plot(p11d2, 'genetic_umap_signatures_per_cluster_tx.pdf', width=12, height=8, plot.if.exists=TRUE)
    }


    # Signature per sample per cluster
    for (cl in 1:mod$G) {
        n <- length(allsigs[class==cl, unique(sample)])
        p11e <- ggplot(allsigs[class==cl], aes(x=mutsig, weight=count, fill=mutsig)) +
            geom_bar() +
            facet_wrap(~interaction(muttype, center, alias), scale='free', ncol=3) +
            theme(aspect.ratio=1/4, axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
            labs(title=paste("Cluster", cl))
            #facet_grid(center+alias~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
            #theme(aspect.ratio=1/4) + labs(title=paste("Cluster", cl))
        do.plot(p11e, paste0('genetic_umap_signatures_per_sample.cluster_', cl, '.pdf'), width=12, height=2+n, plot.if.exists=TRUE)

        if (do.tx) {
            p11e2 <- ggplot(allsigs.tx[class==cl], aes(x=mutsig, weight=count, fill=mutsig)) +
                geom_bar() +
                facet_grid(center+alias~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
                scale_fill_manual(values=tx_cols(c(sbs96_cols(), id83_cols(), dbs78_cols())), guide='none') +
                theme(aspect.ratio=1/4) + labs(title=paste("Cluster", cl))
            do.plot(p11e2, paste0('genetic_umap_signatures_per_sample_tx.cluster_', cl, '.pdf'), width=12, height=2+n, plot.if.exists=TRUE)
        }
    }

    #allsigs[, c('V1', 'V2', 'class') := list(NULL,NULL,NULL)]
    #allsigs.tx[, c('V1', 'V2', 'class') := list(NULL,NULL,NULL)]
}
stop('done')


# Cluster by COSMIC signatures
if (FALSE) {
    # Use COSMIC assignments
    read.sigs <- function(fname) {
        x <- fread(fname)
        x <- x[, .SD, .SDcols=c('Samples',names(x)[-1][colSums(x[,-1])>0])]
        colnames(x)[1] <- 'sample'
        x
    }
    snv.cosmic <- read.sigs('mutsigs/SBS96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
    id.cosmic <- read.sigs('mutsigs/ID83/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
    dbs.cosmic <- read.sigs('mutsigs/DBS78/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')
    all.cosmic <- meta[snv.cosmic[id.cosmic[dbs.cosmic,,on=.(sample)],,on=.(sample)],,on=.(sample)]
    all.cosmic <- all.cosmic[sample != '5871-Neuron-5']
    m <- as.matrix(all.cosmic[, -(1:12)])
    rownames(m) <- all.cosmic[[1]]
    u <- as.data.table(umap(m)$layout, keep.rownames='sample')
    bic <- mclustBIC(u[,.(V1,V2)])
    mod <- Mclust(u[,.(V1,V2)], x=bic)
    u[, class := mod$classification]
    z <- meta[u,,on=.(sample)]
    allsigs <- u[allsigs,,on=.(sample)]
    allsigs.tx <- u[allsigs.tx,,on=.(sample)]

    # "Genetic UMAP"
    p11c <- ggplot(z, aes(x=V1, y=V2, shape=tissue)) +
        geom_mark_hull(aes(group=class, label=class), concavity=10, fill='lightgrey', alpha=0.25) +
        geom_point(aes(col=center), size=3) +
        base.theme + smallerfont
    do.plot(p11c, 'genetic_umap_COSMIC_allsigs.pdf', width=9, height=9)

    # Signatures per cluster
    p11d <- ggplot(allsigs, aes(x=mutsig, weight=w, fill=mutsig)) +
        geom_bar() +
        facet_grid(class~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
        scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
        theme(aspect.ratio=1/4)
    do.plot(p11d, 'genetic_umap_COSMIC_signatures_per_cluster.pdf', width=12, height=8)

    p11d2 <- ggplot(allsigs.tx, aes(x=mutsig, weight=w, fill=mutsig)) +
        geom_bar() +
        facet_grid(class~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
        scale_fill_manual(values=tx_cols(c(sbs96_cols(), id83_cols(), dbs78_cols())), guide='none') +
        theme(aspect.ratio=1/4)
    do.plot(p11d2, 'genetic_umap_COSMIC_signatures_per_cluster_tx.pdf', width=12, height=8)


    # Signature per sample per cluster
    for (cl in 1:mod$G) {
        n <- length(allsigs[class==cl, unique(sample)])
        p11e <- ggplot(allsigs[class==cl], aes(x=mutsig, weight=w, fill=mutsig)) +
            geom_bar() +
            facet_grid(center+alias~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
            scale_fill_manual(values=c(sbs96_cols(), id83_cols(), dbs78_cols()), guide='none') +
            theme(aspect.ratio=1/4) + labs(title=paste("Cluster", cl))
        do.plot(p11e, paste0('genetic_umap_COSMIC_signatures_per_sample.cluster_', cl, '.pdf'), width=12, height=2+n)

        p11e2 <- ggplot(allsigs.tx[class==cl], aes(x=mutsig, weight=w, fill=mutsig)) +
            geom_bar() +
            facet_grid(center+alias~factor(muttype, levels=c('snv', 'indel','dbs')), scale='free') +
            scale_fill_manual(values=tx_cols(c(sbs96_cols(), id83_cols(), dbs78_cols())), guide='none') +
            theme(aspect.ratio=1/4) + labs(title=paste("Cluster", cl))
        do.plot(p11e2, paste0('genetic_umap_COSMIC_signatures_per_sample_tx.cluster_', cl, '.pdf'), width=12, height=2+n)
    }
    allsigs[, c('V1', 'V2', 'class') := list(NULL,NULL,NULL)]
    allsigs.tx[, c('V1', 'V2', 'class') := list(NULL,NULL,NULL)]
}

stop('done')

# PTA can identify through sigs and phylogenies:
# 1. "rare" cell states (i.e., smoking+ cells in an ex-smoker)
# 2. cell types (i.e., SBS1+ colon cells which are likely crypts vs. others)

# shared calls

### STORIES FOR LUNG:
# 1. verifyBAMID = QC
# 2. multi-modal signatures of smoking: SBS, ID, DBS
# 3. what to do with very diff burden estimates? regress on QC to see if it corrects?
# 4. shared mutations - ??? why do they cluster by center

# verifyBAMID2 contamination
p11 <- ggplot(meta[v,,on=.(sample)], aes(x=fct_reorder(sample, .x=FREEMIX), y=FREEMIX, fill=tissue)) +
    geom_col() +
    facet_grid(~center+tissue, scale='free_x', space='free_x') +
    theme(axis.text.x=element_text(angle=90, vjust=1/2, hjust=1)) +
    xlab(element_blank()) + ylab("Contamination ('FREEMIX')")
do.plot(p11, 'verifybamid2.pdf', width=16, height=6)

# MAPD explains most of verifybamid2 contamination - a few outliers
z <- melt.and.annotate(mapd(sm))[binsize==64000][, mapd := value][,binsize := NULL][,value:=NULL][]
p11b <- ggplot(z, aes(x=mapd, y=FREEMIX, col=center, label=alias)) +
    geom_point() +
    geom_smooth(data=z[!(mapd>0.6 |FREEMIX>0.2)], aes(group=1), se=F, col='black') +
    geom_text_repel(data=z[mapd>0.6 |FREEMIX>0.2]) +
    ylab('Contamination ("FREEMIX")') + xlab("MAPD") +
    base.theme + bigfont + theme(aspect.ratio=1)
do.plot(p11b, 'verifybamid2_vs_mapd.pdf', width=7, height=7)

if (FALSE) {
    # Try beta fitting to summarize AB distribution in one number
    # This takes a while (~5-10 min).
    # Directionally correct - higher shape1 and shape2 = more well-balanced
    # but the fit is not particularly good because of the masses at 0 and 1.
    # Need a mixture.
    # Is very closely related to ADO = mean(af <= 0.01 | af >= 0.99)
    system.time(bfits <- rbindlist(lapply(sm, function(s) {
        print(name(s))
        # pmin(pmax()) is critical: dbeta(0) and dbeta(1) = 0
        af <- training(s)[muttype=='snv' & bulk.gt=='0/1'  & bulk.dp>64 & dp>15, pmin(1-1e-4, pmax(1e-4, af))]
        f <- MASS::fitdistr(x=af, densfun='beta', start=list(shape1=1/2,shape2=1/2))
        data.table(sample=name(s), t(f$estimate))
    })))
    meta <- meta[bfits[,.(sample,ab.beta=shape1)],,on=.(sample)]
}
# IF TIME:
# replace N calls with digested N calls
# add plot showing number of calls filtered by digest
# add plot showing number of added rescue calls

# chr4 on SM-OGQSQ example: ~0 depth chrom except small region, in small region clustered somatic
# calls that are normally removed in a post-processing "clustered" filter.
