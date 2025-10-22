#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 240:00:00
#SBATCH -p priopark
#SBATCH -A park_contrib

outdir=panel_104_PTA_2_bulk_ST002_52_PTA_18_bulk_neurons

if [ 0 -eq 1 ]; then
/n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/bin/scan2 -d $outdir \
    --snakefile /n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/snakemake/Snakefile init

    # Old, non-DAC ref
    #--ref /n/data1/hms/dbmi/park/jluquette/scan2-test/hg38/Homo_sapiens_assembly38.fasta \
    # NOTE: most of the makepanel run used the wrong, above reference. However, the primary
    # contigs are identical to the DAC reference, so there should be no difference here
    # (unlike in alignment, where presence of other contigs can change where reads align).
/n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/bin/scan2 -d $outdir \
    --snakefile /n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/snakemake/Snakefile config \
    --verbose \
    --analysis makepanel \
    --gatk sentieon_joint \
    --genome hg38 \
    --ref /n/data1/hms/dbmi/park-smaht_dac/ref/GRCh38_no_alt/hg38_no_alt.fa \
    --dbsnp /n/data1/hms/dbmi/park/jluquette/scan2-test/hg38/common_all_20180418.chrprefix.vcf \
    --scripts /n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/scripts \
    --resources /n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/resources \
    --regions-file /n/data1/hms/dbmi/park/jluquette/SMaHT/pta_pilot/scan2_analysis_regions_no_PARs_chr1-22XY_30839windows_100kb.txt \
    --makepanel-metadata /n/data1/hms/dbmi/park/jluquette/SMaHT/pta_pilot/panel_metadata.csv \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_bulkWgs/seq_data/SMAFIGBTW73V.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_bulkWgs/seq_data/SMAFI43OG9OX.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI1HFQVKJ_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI1U6CPYI_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI22SJWPL_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI2OGACJE_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI2S6UC16_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI2YMLDA7_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI31YVJ3T_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI33GSPI6_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI3AU2533_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI3D883AK.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI3LQMZAM.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI3MYQUKI_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI5PYYVZA_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI68KX6AU_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI7B7LEKT_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI7DSNDAI_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI7IMXNTK_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI7LUMPMI_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI7TOOAK6_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI8665PKU_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI8INXOJ2.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI8RHUJVS_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFI8TN14MW_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI8XGLKA8.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI8ZRW634_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI97SY757_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI981BRET_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI9KTPSTS_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFI9VX7TP6_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIAL5DGVG.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIBGR661P.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIBRG71BL_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIC2GZU55_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIC4JIWXS_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFICVRG2S2.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFID1IBZZW_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIDDIVNII.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIDHFGEO7.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIEEVQFDA_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIEP1ZEN3.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIER1L2SB_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIER63KAS_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIFNLHOX2_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIFUJE1HS_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIG9BD627_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIGL3ODC4_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIGS48K1K_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIH7BK4MN_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIHS8OJ7V_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIIC9KFU3_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIJ1NNS51_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIJ38GDVW_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIJ5TQYKQ_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIJ99U4UM_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIJB1EPU9_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIJNLWWF8_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIK2DTD5D_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIKDLL4JW.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIKHXJF8Y_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIKJ42GM9_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIKWGOS4G_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIL87OPFP_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIMBNSZQJ.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIMHL5HPA.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIMKBOQDI_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIMYRD1PF_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIMZPF1D1_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIN8PYW1B_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIN8UINE3_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFING6ISUC_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIO2UH81W_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIOGMPQW6_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIOJS918Q_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIOOY4YCE_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIP9OVY1K_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIPAETAVA_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIPFCBZT2_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIPPXLCUW_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIQ9M2ZN3_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIQEEU8KJ_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIRITNLAP_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFITFI5XDC_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFITKGMLQC_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFITNKFQ7I_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFITZ6QF7U_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIURQ1S35_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIUXVN8X8_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIV5JQ1RF_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIVPONF1W_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIVSA5AI5_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_BCM/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIWA7N6CM_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIWCKZ4KG_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIWPEW7YV_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_NYGC/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIWZVXVSI.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIXH65Z9X_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIXO7Q97U_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIXX3VOOS_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIZ8WLLK3_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/TTD_Mayo/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIZBCMB9O.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1G/illuminaNovaseq_pta/seq_data/SMAFIZLPQFUB_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_Broad/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIZMGF89R_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park-smaht_dac/DATA/GCC_WashU/ST002-1D/illuminaNovaseq_pta/seq_data/SMAFIZYIW8RA_tags_removed.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1278BA9-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1278BA9-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1278BA9-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1278_heart_bulk_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1465-cortex_BulkDNA_WGSb_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1465-heart_BulkDNA_WGSb_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1465BA9-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1465BA9-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1465BA9-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/1465BA9-D_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4638-Bulk-Heart_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4638-Neuron-4_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4638-Neuron-5_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4638-Neuron-6_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4643-Neuron-3_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4643-Neuron-4_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4643-Neuron-6_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/4643_Bulk-Liver_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5087-hrt-1b1-final-all_realigned_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5087PFC-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5087PFC-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5087PFC-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5219-Neuron-2_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5219-Neuron-4_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5219-Neuron-5_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5219_cb_bulk_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5559PFC-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5559PFC-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5559PFC-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5559_0928-hrt-1b1_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5657PFC-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5657PFC-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5657PFC-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5657_0717-hrt-1b1_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5817PFC-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5817PFC-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5817PFC-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5817_liver_bulk_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5823-tempmusc-1b1_20170221-WGS-final-all_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5823PFC-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5823PFC-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5823PFC-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5871-BLK-liver_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5871-Neuron-4_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5871-Neuron-5_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/5871-Neuron-6_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/936-hrt-1b1_20170221-WGS-final-all_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/936PFC-A_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/936PFC-B_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/936PFC-C_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB4976_E1_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB4976_E2_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB4976_E3_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB4976_bulk_PTA_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5451_B2_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5451_B3_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5451_B5_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5451_bulk_PTA_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5572_D2_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5572_D3_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5572_D4_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5572_bulk_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5666_F1_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5666_F2_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5666_F5_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5666_bulk_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5943_C2_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5943_C4_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5943_C5_qualcal.bam \
    --bam /n/data1/hms/dbmi/park/jbrew/joe_realignment_no_decoy/recal/UMB5943_bulk_PTA_qualcal.bam


/n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/bin/scan2 -d $outdir \
    --snakefile /n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/snakemake/Snakefile \
    validate
fi

/n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/bin/scan2 -d $outdir  \
    --snakefile /n/data1/hms/dbmi/park/jluquette/scan2-test/SCAN2/snakemake/Snakefile \
    makepanel \
        --joblimit 5000 \
        --n-cores 8 \
        --snakemake-args ' --set-resources gatk_scatter:mem_mb=100000 --notemp --keep-going --latency-wait=60 --rerun-incomplete --executor slurm --default-resources slurm_account=park slurm_partition=short,park runtime=720 tmpdir="/n/data1/hms/dbmi/park/jluquette/SMaHT/pta_pilot/tmp" --max-jobs-per-timespan=6/1s --max-status-checks-per-second=0.1'
