#plot.binned.counts(sm[c(cm[center=='mayo-washu' & alias=='Lu_14', sample], cm[center=='g1' & alias=='Lu-12',sample], cm[center=='bch-broad' & alias=='Co_9',sample], cm[center=='mayo-washu' & alias=='Lu_n4',sample])], type='count')

four.cells <- c(cm[center=='mayo-washu' & alias=='Lu_14', sample], cm[center=='g1' & alias=='Lu-12',sample], cm[center=='bch-broad' & alias=='Co_9',sample], cm[center=='mayo-washu' & alias=='Lu_n4',sample])

pdf('supp_figure10.pdf', width=5.0, height=5.75)
layout(1:4)
for (cn in four.cells)
    helper.plot.binned.counts(collapse.binned.counts(sm[[cn]]@binned.counts$sc, 5), type='count', sample.name=cn)
dev.off()
