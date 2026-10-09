# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

dmr_path=file.path(UPDATE_ROOT, "整合数据", "0. MP_4215-SZ-DMR.txt")
dmr=read.table(dmr_path,head=T,sep="\t")

gwas_loci_path=file.path(UPDATE_ROOT, "整合数据", "PGC3_SZ_Loci.txt")
gwas_loci=read.table(gwas_loci_path,head=T,sep="\t")

countall=c()
for (i in 1:nrow(dmr)){
  chr=gsub("chr","",as.character(dmr[i,1]))
  st=dmr[i,2]
  ed=dmr[i,3]
  count=which(gwas_loci[,1]==chr & gwas_loci[,11]>=st & gwas_loci[,10]<=ed)
  if (length(count)==0){
    count=c(count,NA)
  }
  if (length(count)==1){
    countall=c(countall,count)
  }
  if (length(count)>1){
    print("error")
  }
}

final=cbind(dmr,gwas_loci[countall,])

write.table(final,file.path(UPDATE_ROOT, "整合数据", "1.DMR进行GWAS_Loci注释.txt"),quote=F,row.names=F,sep="\t")

