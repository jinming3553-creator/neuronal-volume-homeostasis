# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

dmr_path=file.path(UPDATE_ROOT, "结果", "1.DMR与GWAS交集部分.txt")
dmr=read.table(dmr_path,head=T,sep="\t")

probe_path=file.path(UPDATE_ROOT, "结果", "2.PGC3-Primary-mQTL-SMR合并DMP.txt")
probe=read.table(probe_path,head=T,sep="\t")


count1=c()
count2=c()
for (i in 1:nrow(dmr)){
  chr=gsub("chr","",as.character(dmr[i,1]))
  st=dmr[i,2]-2500
  ed=dmr[i,3]+2500
  count=which(probe[,2]==chr & probe[,4]>=st & probe[,4]<=ed)
  count1=c(count1,rep(i,length(count)))
  count2=c(count2,count)
}

final=cbind(dmr[count1,],probe[count2,])

write.table(final,file.path(UPDATE_ROOT, "结果", "3.DMR_in_PGC3_Loci与Probe交集.测试.txt"),quote=F,row.names=F,sep="\t")


