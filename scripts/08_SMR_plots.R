# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

source(file.path(EXTERNAL_CODE_ROOT, "plot_SMR.r")) 

# Read the data file in R:
SMRData = ReadSMRData(file.path(UPDATE_ROOT, "结果", "LRRC8A.cg13630871.mMeta.SMR作图.txt"))
# Plot the SMR results in a genomic region centred around a probe:
SMRLocusPlot(data=SMRData, smr_thresh=1E-4, heidi_thresh=0.05, plotWindow=200)
# smr_thresh: genome-wide significance level for the SMR test.
# heidi_thresh: threshold for the HEIDI test. The default value is 0.05.
# cis_wind: size of a window centred around the probe to select cis-eQTLs for plot. The default value is 2000Kb.
# max_anno_probe: maximum number of probe names to be displayed on the figure. The default value is 16.
# Plot effect sizes from GWAS against those from eQTL study:
SMREffectPlot(data=SMRData, trait_name="SCZ") 
# trait_name: name of the trait or disease. 



SMRData = ReadSMRData(file.path(UPDATE_ROOT, "结果", "LRRC8A_eMeta.ENSG00000136802.12作图.txt"))
SMRLocusPlot(data=SMRData, smr_thresh=1E-4, heidi_thresh=0.05, plotWindow=200)


