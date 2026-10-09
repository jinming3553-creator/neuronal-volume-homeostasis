# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

load(file.path(DATA_ROOT, "5. CpG-Gene关联", "allgene.save"))
load(file.path(DATA_ROOT, "5. CpG-Gene关联", "BrainSeq_PhaseI_Gene.Fpkm.save"))

loc=file2[which(file2[,1]=="9" & file2[,2]>=130583092 & file2[,2]<=132583092),]
id=intersect(rownames(fpkm),as.character(loc[,4]))

rownames(loc)=loc[,4]
loc=loc[id,]
genes=t(fpkm[id,])

write.table(loc,file.path(DATA_ROOT, "5. CpG-Gene关联", "Annotation_for_genes_around.txt"),quote=F,row.names=F,sep="\t")
write.table(genes,file.path(DATA_ROOT, "5. CpG-Gene关联", "BrainSeq_PhaseI_Gene.around.Fpkm.txt"),quote=F,row.names=F,sep="\t")

pheno_geo=read.table(file.path(DATA_ROOT, "5. CpG-Gene关联", "phenotype_data_forGeo.csv"),head=T,sep=",")
pheno_dabap=read.table(file.path(DATA_ROOT, "5. CpG-Gene关联", "phenotype_annotated_szControlEqtl_DLPFC_withDBGAP.csv"),head=T,sep=",")
pheno=merge(pheno_geo[,1:2],pheno_dabap[,c(1,2,5)],by="BrNum")
cpg=read.table(file.path(DATA_ROOT, "5. CpG-Gene关联", "GSE74193_GEO_procData.cg13630871.csv"),head=T,sep=",")
colnames(cpg)[1]="SampleID"
pheno_cpg=merge(pheno,cpg[,c(1,3,2)],by="SampleID")
#rownames(pheno_cpg)=pheno_cpg[,3]
pheno_cpg=pheno_cpg[which(pheno_cpg$Age>0),]
id=unique(as.character(pheno_cpg[,3]))

result=matrix(,nrow=length(id),ncol=6)
for (i in 1:length(id)){
  rms=mean(pheno_cpg[which(pheno_cpg[,3]==id[i]),6])
  result[i,1:5]=as.character(pheno_cpg[which(pheno_cpg[,3]==id[i])[1],1:5])
  result[i,6]=rms
}
colnames(result)=colnames(pheno_cpg)
pheno_cpg=data.frame(result)
rownames(pheno_cpg)=pheno_cpg[,3]

RNum=intersect(rownames(genes),as.character(pheno_cpg[,3]))
cpg_genes=cbind(pheno_cpg[RNum,],genes[RNum,])


result=matrix(,nrow=109,ncol=3)
rms=as.numeric(cpg_genes[,6])
for (i in 1:109){
  fpkm=cpg_genes[,i+6]
  result[i,3]=cor.test(rms,fpkm)$p.value
  result[i,2]=cor(rms,fpkm)
}
result[,1]=colnames(cpg_genes)[7:ncol(cpg_genes)]
colnames(result)=c("TargetID","r","P")

final=merge(loc,result,by="TargetID")
final=final[order(final$start),]



