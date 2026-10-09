# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

dmr_path=file.path(UPDATE_ROOT, "0. Twin-DMR", "0. MP_4215-SZ-DMR.txt")
dmr=read.table(dmr_path,head=T,sep="\t")

mqtl_epi_path=file.path(UPDATE_ROOT, "4. meQTL-GWAS SMR", "Brain-mMeta.epi")
meqtl_epi=read.table(mqtl_epi_path,head=F,sep="\t")

final=c()
for (i in 1:22){
  selfile2=meqtl_epi[which(meqtl_epi[,1]==i),]
  selfile1=dmr[which(dmr[,1]==paste("chr",i,sep="")),]
  count=c()
  for (j in 1:nrow(selfile1)){
    st=selfile1[j,2]-2500
    ed=selfile1[j,3]+2500
    count=c(count,which(selfile2[,4]<=ed & selfile2[,4]>=st))
  }
  count=unique(count);count=count[order(count)]
  final=rbind(final,selfile2[count,])
}
write.table(final,file.path(UPDATE_ROOT, "4. meQTL-GWAS SMR", "Brain-mMeta.DMR.epi"),quote=F,row.names=F,sep="\t")
