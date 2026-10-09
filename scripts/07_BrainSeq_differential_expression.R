# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

final=read.table(file.path(DATA_ROOT, "6. DEG", "BrainSeq_PhaseII_LRRC8A表达数据.txt"),head=T,sep="\t")

file1=final[which(final$Region=="DLPFC"),]
file2=final[which(final$Region=="HIPPO"),]

x1=file1$Dx
x2=file2$Dx


result=matrix(,nrow=5,ncol=7)
for (i in 9:ncol(file1)){
  y1=file1[,i]
  y2=file2[,i]
  result[i-8,2]=median(y1[which(x1=="Control")])
  result[i-8,3]=median(y1[which(x1=="Schizo")])
  result[i-8,4]=wilcox.test(y1[which(x1=="Control")],y1[which(x1=="Schizo")])$p.value
  result[i-8,5]=median(y2[which(x2=="Control")])
  result[i-8,6]=median(y2[which(x2=="Schizo")])
  result[i-8,7]=wilcox.test(y2[which(x2=="Control")],y2[which(x2=="Schizo")])$p.value
}
result=data.frame(result)
result[,1]=colnames(file1)[9:13]

