# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

dmr = read.table(file.path(UPDATE_ROOT, "结果", "1.DMR进行GWAS_Loci注释.txt"),head=T,sep="\t")
dmp = read.table(file.path(UPDATE_ROOT, "结果", "2.PGC3-Primary-mQTL-SMR合并DMP.txt"),head=T,sep="\t")

dmr = dmr[which(abs(dmr[,11])>1 & dmr[,15]<0.01),]
chr=gsub("chr","",as.character(dmr[,1]))
st=dmr[,2]-2500
ed=dmr[,3]+2500

count=c()
for (i in 1:length(chr)){
  count=c(count,which(dmp[,2]==chr[i] & dmp[,4]>=st[i] & dmp[,4]<=ed[i]))
}

count=unique(count)
count=count[order(count)]

dmp=dmp[count,]
write.table(dmp,file.path(UPDATE_ROOT, "结果", "2.PGC3-Primary-mQTL-SMR合并DMP与显著DMR交集.txt"),quote=F,row.names=F,sep="\t")
