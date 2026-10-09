# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

dmp_path=file.path(UPDATE_ROOT, "3. DMP", "jaffe_onlineTable1_DMPsWithExprs.csv")
dmp=read.table(dmp_path,head=T,sep=",")

mqtl_smr_path=file.path(UPDATE_ROOT, "4. meQTL-GWAS SMR", "PGC3_Primary_GWAS_Brain_mMeta_meQTL_SMR.txt")
mqtl_smr=read.table(mqtl_smr_path,head=T,sep="\t")


rownames(dmp)=dmp[,1]
rownames(mqtl_smr)=mqtl_smr[,1]
id=intersect(rownames(dmp),rownames(mqtl_smr))

final=cbind(mqtl_smr[id,],dmp[id,2:6])

write.table(final,file.path(UPDATE_ROOT, "结果", "2.PGC3-Primary-mQTL-SMR合并DMP.txt"),quote=F,row.names=F,sep="\t")


