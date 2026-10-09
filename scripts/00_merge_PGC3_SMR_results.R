# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

pri_path=file.path(UPDATE_ROOT, "整合数据", "PGC3_Primary.txt")
cor_path=file.path(UPDATE_ROOT, "整合数据", "PGC3_CORE.txt")
eur_path=file.path(UPDATE_ROOT, "整合数据", "PGC3_EUR.txt")


pri=read.table(pri_path, header=TRUE, sep="\t")
core=read.table(cor_path, header = TRUE, sep="\t")
eur=read.table(eur_path,header = TRUE, sep="\t")

rownames(pri)=pri[,1]
rownames(eur)=eur[,1]
rownames(core)=core[,1]

length(intersect(rownames(core),rownames(pri)))
length(intersect(rownames(core),rownames(eur)))
length(intersect(rownames(eur),rownames(pri)))

id=unique(c(rownames(core),rownames(eur),rownames(pri)))

mqtlsmr=cbind(pri[id,c("probeID","p_SMR","p_SMR_multi")],core[id,c("p_SMR","p_SMR_multi")],eur[id,c("p_SMR","p_SMR_multi")])
mqtlsmr[,1]=id
colnames(mqtlsmr)=c("probeID","Pri_P_SMR","Pri_P_SMR_Multi","Core_P_SMR","Core_P_SMR_Multi","EUR_P_SMR","EUR_P_SMU_Multi")

write.table(mqtlsmr,file.path(UPDATE_ROOT, "整合数据", "PGC3_mQTL_SMR.txt"),quote=F,row.names=F,sep="\t")


