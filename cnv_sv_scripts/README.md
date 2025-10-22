# CNV/SV Scripts

This repository contains utilities for calling, genotyping, and visualising structural variants (SVs) and copy-number variation (CNV) signals across single-cell PTA datasets. Below is a quick tour of each script along with dependency notes.

## Requirements

- Python 3.9+
- Bash shell
- [CNVpytor v1.3.2](https://github.com/abyzovlab/CNVpytor/releases/tag/v1.3.2)
- [Manta v1.6.0](https://github.com/Illumina/manta/releases/tag/v1.6.0)
- Python libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`, `cnvpytor`, `matplotlib`, `nbformat`

Ensure CNVpytor `.pytor` files and Manta VCF outputs are organised as expected by the scripts (see individual sections).

## Scripts

### `run_cnvpytor.sh`
Wrapper around CNVpytor that builds histograms, B-allele frequencies, and CNV calls from a BAM and phased VCF. Accepts an optional core count (defaults to 8) and produces `.pytor` files plus TSV call files in the working directory.

Before running:
- Provide aligned single cell BAM and phased bulk VCF paths as the first two arguments.
- Ensure CNVpytor v1.3.2 is installed.
- Adjust the optional cores argument if you want something other than the default of 8.

### `run_manta.sh`
Configures and launches the Manta somatic workflow for paired single cell/blood BAMs. The script validates inputs, creates the run directory, and runs the workflow with a configurable core count.

Before running:
- Supply single cell PTA BAM, bulk BAM, reference FASTA, and a run directory.
- Confirm Manta v1.6.0 is installed.
- Optionally override the default core count (8) to match your hardware.

### `filtering_cnv_calls.py`
Converts CNVpytor call tables (`calls/*.2d.100k.tsv`) into an aggregated DataFrame, tags large events (deletion, duplication, CNN-LOH), and exports summary plots and statistics (`figures/`). Also writes a `large_calls.tsv` table with merged segments.

Before running:
- Populate `calls/` with CNVpytor call tables (e.g. `sample.2d.100k.tsv`).
- Provide a `samples_tissue` mapping file containing sample IDs and tissue codes ("C" for colon and "L" for lung)
- Create (or let the script create) a writable `figures/` directory for plots.

### `genotyping_manta_calls.py`
Parses Manta somatic VCFs, summarises DEL/DUP counts per sample, and genotypes events against single-cell CNVpytor data to estimate copy number (`rd`). Outputs histograms, violin plots, and filtered call tables.

Before running:
- Ensure `samples_tissue` lists every sample to analyse with its tissue type.
- Place Manta output files in `manta_vcfs/` as `<sample>.somaticSV.vcf.gz`.
- Store matching CNVpytor roots in `pytor/<sample>.pytor` for rd genotyping.
- Verify CNVpytor v1.3.2 is installed; without it, rd values fall back to empty.

### `t_cells.py`
Focuses on T-cell receptor loci: aggregates Manta SVs, computes CNVpytor read-depth tracks, fits discrete CN states, and renders step and heatmap panels for control vs. T-cell samples (`figures/figC_*.png`).

Before running:
- Edit `CONTROL_SAMPLE` and `T_CELLS` inside the script to your actual sample IDs.
- Provide matching Manta VCFs (`manta_vcfs/`) and CNVpytor roots (`pytor_phased/`).
- Confirm `samples_tissue` includes all samples referenced in the constants.

### `phase_chrY.py`
Phases SNPs within the PAR1 region using CNVpytor outputs from samples with chromosome Y loss. Generates coverage histograms, heatmaps, and BAF distributions for loss/gain/control cohorts.

Before running:
- Verify `ids` and `ids_loss` arrays match the sample IDs you wish to analyse.
- Place corresponding CNVpytor root files in `pytor/<sample>.pytor`.
- Make sure CNVpytor v1.3.2 and matplotlib are available in your Python environment.
