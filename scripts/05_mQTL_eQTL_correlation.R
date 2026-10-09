# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

eqtl=read.table(file.path(DATA_ROOT, "5. mQTL-eQTL", "LRRC8A_CpG附近基因.eqtl.txt"),head=T,sep="\t")
mqtl=read.table(file.path(DATA_ROOT, "5. mQTL-eQTL", "LRRCA_CpG.mqtl.txt"),head=T,sep="\t")
eqtl=eqtl[which(eqtl[,2]=="9"),]

cpg=mqtl[which(mqtl$featureName=="cg13630871"),]

gid=unique(as.character(eqtl$featureName))

result=matrix(,nrow=length(gid),ncol=3)
for (i in 1:length(gid)){
  target=eqtl[which(eqtl[,4]==gid[i]),]
  final=merge(cpg,target,by="SNPid")
  #if (nrow(final)>1){
  #final=final[which(final$pValue.x<0.05),]
  if (nrow(final)>=10){
  x=final$SpearmanRho.x
  y=final$SpearmanRho.y
  result[i,3]=cor.test(x,y)$p.value
  result[i,2]=cor(x,y)
  #}
  }
}
result=data.frame(result)
result[,1]=gid
result=result[order(result[,3]),]
colnames(result)=c("Gene.name","r","P")

annotation=read.table(file.path(DATA_ROOT, "5. mQTL-eQTL", "LRRC8A_CpG附近基因.txt"),head=T,sep="\t")
final=merge(annotation[,c(1,2,10,5,6)],result,by="Gene.name")
write.table(final,file.path(DATA_ROOT, "5. mQTL-eQTL", "LRRC8A_CpG附近基因-mQTL-eQTL相关.txt"),quote=F,row.names=F,sep="\t")

final=read.table(file.path(DATA_ROOT, "结果", "LRRC8A_CpG附近基因-mQTL-eQTL相关.txt"),head=T,sep="\t")
x=final$TSS
y=-log10(final$P)

plot(x,y,pch=16,xlab="Position",ylab="-log10(Pvalue)",main="Correlation between effect of mQTL and eQTL from Brain-xQTL dataset")
abline(h=100,col="red",lwd=2,lty=2)

final=final[order(final$P),]
gid=as.character(final[which(final$P<1E-100),1])

for (i in 1:length(gid)){
  target=eqtl[which(eqtl[,4]==gid[i]),]
  final=merge(cpg,target,by="SNPid")
  #if (nrow(final)>1){
  #final=final[which(final$pValue.x<0.05),]
  if (nrow(final)>=10){
    x=final$SpearmanRho.x
    y=final$SpearmanRho.y
    plot(x,y,pch=16,xlab="Spearman Rho of mQTL",ylab="Spearman Rho of eQTL",main=gid[i])
    abline(lm(y~x),col="red",lwd=2,lty=2)
    #}
  }
}
result=data.frame(result)
result[,1]=gid
result=result[order(result[,3]),]
colnames(result)=c("Gene.name","r","P")



