# Script to phase SNPs within the PAR1 region.
# It uses samples exhibiting chromosome Y loss to derive phasing,
# evaluates phasing agreement across samples with loss, and visualizes allele
# frequency distributions for loss, gain, and control samples.
#
# Requirements: CNVpytor pytor files in `pytor/`, `cnvpytor`, `matplotlib`,
# and `numpy` Python packages.


import cnvpytor
import matplotlib.pyplot as plt
import numpy as np


ids_loss = [
    "SMAFIKHXJF8Y",
    "SMAFI22SJWPL",
    "SMAFIFNLHOX2",
    "SMAFIAL5DGVG",
    "SMAFI7B7LEKT",
    "SMAFI9KTPSTS",
    "SMAFIVPONF1W",
]

ids = ids_loss + [
    "SMAFI5PYYVZA", #chrY gain
    "SMAFI7IMXNTK", #chrY gain
    "SMAFIFUJE1HS"  #control cell
]

view = cnvpytor.Viewer([f"pytor/{id}.pytor" for id in ids_loss], {"bin_size": 100000})

snp_alt_zero = {}
snp_ref_zero = {}
cov = []
snp_nonP_alt_zero = {}
snp_nonP_ref_zero = {}
cov_nonP= []

for i in range(len(ids_loss)):
    pos, ref, alt, nref, nalt, gt, flag, qual = view.io[i].read_snp("chrX")
    for p, r, a, n_r, n_a, g, f, q in zip(pos, ref, alt, nref, nalt, gt, flag, qual):
        if ((g%4) in {1,2}) and p<2800000:
            if f == 0:
                cov_nonP.append(n_a+n_r)
                if n_a==0 and n_r > 0:
                    if (p, r, a) not in snp_nonP_alt_zero:
                        snp_nonP_alt_zero[(p,r,a)] = []
                    snp_nonP_alt_zero[(p,r,a)].append(i)
                if n_r==0 and n_a > 0:
                    if (p, r, a) not in snp_nonP_ref_zero:
                        snp_nonP_ref_zero[(p,r,a)] = []
                    snp_nonP_ref_zero[(p,r,a)].append(i)
            else:    
                cov.append(n_a+n_r)
                if n_a==0 and n_r > 0:
                    if (p, r, a) not in snp_alt_zero:
                        snp_alt_zero[(p,r,a)] = []
                    snp_alt_zero[(p,r,a)].append(i)
                if n_r==0 and n_a > 0:
                    if (p, r, a) not in snp_ref_zero:
                        snp_ref_zero[(p,r,a)] = []
                    snp_ref_zero[(p,r,a)].append(i)



plt.rcParams["figure.figsize"] = (8, 10)
plt.hist([cov, cov_nonP], bins=range(0, 51), align='left', rwidth=0.8, stacked=True, label=['P region', 'non P region'])
plt.xlabel("Coverage")
plt.ylabel("Number of HETs")
plt.legend()
plt.show()

plt.rcParams["figure.figsize"] = (4, 8)
hm_a = []
hm_nonP_a = []
for (p, r, a) in snp_nonP_alt_zero:
    hm_nonP_a.append(len(snp_nonP_alt_zero[(p, r, a)]))
for (p, r, a) in snp_alt_zero:
    hm_a.append(len(snp_alt_zero[(p, r, a)]))
plt.hist([hm_a, hm_nonP_a], bins=range(1, 9), orientation='horizontal', align='left', rwidth=0.8, stacked=True, label=['P region', 'non-P region'])
plt.ylabel("#samples sharing zero-ALT-count HET", fontsize=20)
plt.xlabel("#HETs", fontsize=20)
#bottom right corner
plt.legend(loc='lower right',fontsize=16)
plt.xticks([0,200,400,600,800],fontsize=16)
plt.yticks(fontsize=16)
plt.grid()
plt.show()

plt.rcParams["figure.figsize"] = (4, 8)
hm_r = []
hm_nonP_r = []
for (p, r, a) in snp_nonP_ref_zero:
    hm_nonP_r.append(len(snp_nonP_ref_zero[(p, r, a)]))
for (p, r, a) in snp_ref_zero:
    hm_r.append(len(snp_ref_zero[(p, r, a)]))
plt.hist([hm_r, hm_nonP_r], bins=range(1, 9), orientation='horizontal', align='left', rwidth=0.8, stacked=True, label=['P region', 'non-P region'])
plt.ylabel("#samples sharing zero-REF-count HET", fontsize=20)
plt.xlabel("#HETs", fontsize=20)
plt.legend(loc='lower right', fontsize=16)
plt.xticks([0,200,400,600],fontsize=16)
plt.yticks(fontsize=16)
plt.grid()
plt.show()

phased_snps = {}
c1, c2 = 0, 0
for (p, r, a) in snp_alt_zero:
    if len(snp_alt_zero[(p, r, a)]) > 4:
        phased_snps[(p, r, a)] = 1
        c1 += 1
for (p, r, a) in snp_ref_zero:
    if len(snp_ref_zero[(p, r, a)]) > 4:
        phased_snps[(p, r, a)] = 2
        c2 += 1
for (p, r, a) in snp_nonP_alt_zero:
    if len(snp_nonP_alt_zero[(p, r, a)]) > 4:
        phased_snps[(p, r, a)] = 1
        c1 += 1
for (p, r, a) in snp_nonP_ref_zero:
    if len(snp_nonP_ref_zero[(p, r, a)]) > 4:
        phased_snps[(p, r, a)] = 2
        c2 += 1

print(f"Phased {c1} SNPs with alt=0 and ref>0")
print(f"Phased {c2} SNPs with ref=0 and alt>0")

sbafs = []

for sample in [f"pytor/{id}.pytor" for id in ids]:
    view = cnvpytor.Viewer([sample], {"bin_size": 100000})
    pos, ref, alt, nref, nalt, gt, flag, qual = view.io[0].read_snp("chrX")
    spos = []
    sbaf = []
    sc1 = []
    sc2 = []
    bin_size = 100000
    bins = 2800000 // bin_size
    binned_c1 = np.zeros(int(bins))
    binned_c2 = np.zeros(int(bins))
    for p, r, a, n_r, n_a, g, f, q in zip(pos, ref, alt, nref, nalt, gt, flag, qual):
        if ((g%4) in {1,2}) and p<2800000 and f==2:
            if (p, r, a) in phased_snps:
                spos.append(p/1e6)
                if phased_snps[(p, r, a)] == 1:
                    sc1.append(n_a)
                    sc2.append(n_r)
                    sbaf.append(n_a / (n_a + n_r))
                    binned_c1[int(p/bin_size)] += n_a
                    binned_c2[int(p/bin_size)] += n_r
                else:
                    sc1.append(n_r)
                    sc2.append(n_a)
                    sbaf.append(n_r / (n_a + n_r))
                    binned_c1[int(p/bin_size)] += n_r
                    binned_c2[int(p/bin_size)] += n_a
    sbafs.append(sbaf)

plt.rcParams["figure.figsize"] = (5, 16)
fig, axes = plt.subplots(len(sbafs), 1, sharex=True, figsize=(5, 16))
for i in range(len(sbafs)):
    ax = axes[i] if len(sbafs) > 1 else axes
    ax.hist(sbafs[i], bins=21, range=(0, 1))
    ax.grid(True)
    # Set x label and ticks on top
    if i == 0:
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.set_xticks([0, 1/3, 1/2, 2/3, 1])
        ax.set_xticklabels(['0', '1/3', '1/2', '2/3', '1'], fontsize=12)
    else:
        ax.set_xticklabels([])
    ax.tick_params(axis='y', labelsize=14, right=True, labelright=True, left=False, labelleft=False)
plt.tight_layout(h_pad=0.2)
plt.show()
