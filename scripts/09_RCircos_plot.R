# Reproducibility script for the LRRC8A study
# Cleaned for public release: local absolute paths were replaced by configurable roots.
# Set these environment variables if your data are stored elsewhere.
DATA_ROOT <- Sys.getenv("LRRC8A_DATA_ROOT", unset = "../data/LRRC8A")
UPDATE_ROOT <- Sys.getenv("LRRC8A_UPDATE_ROOT", unset = "../data/LRRC8A_update")
EXTERNAL_CODE_ROOT <- Sys.getenv("LRRC8A_EXTERNAL_CODE_ROOT", unset = "../external_code")

library(ggplot2)
library(RCircos)  # 加载RCircos包

##导入各种需要的数据，感兴趣可以自己打开看看数据是什么样的
data (UCSC.HG19.Human.CytoBandIdeogram)
data(RCircos.Scatter.Data)
data(RCircos.Histogram.Data)
data(RCircos.Link.Data)
data(RCircos.Ribbon.Data)

# 导入内建人类染色体数据
data(UCSC.HG19.Human.CytoBandIdeogram);
#设置不显示的染色体，如 c(1,3)   
chr.exclude <- c("chrY"); 
#设置染色体数据
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram; 
#设置内部环形个数
tracks.inside <- 3; 
#设置外部环形个数
tracks.outside <- 0; 
# 导入上面四个基本参数
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
# 列出所有绘图参数
RCircos.List.Plot.Parameters()
# 定义要生成的图像文件
#pdf(file="RCircosDemoHumanGenome", height=8, width=8, compress=TRUE);
RCircos.Set.Plot.Area();
# height和width指定生成图片的长和宽，compress指定生成的图片是否需要压缩
# 绘制染色体图形
RCircos.Chromosome.Ideogram.Plot();


#####散点图部分-SMR
df = read.table(file.path(UPDATE_ROOT, "结果", "1. PGC3_Primary_GWAS_Brain_mMeta_meQTL_SMR.txt"),head=T,sep="\t")
df$Chrom=paste("chr",df$ProbeChr,sep="")
df$Start=df$Probe_bp
df$Stop=df$Start+1
df$logSMR=-log10(df$p_SMR)
smr=df[,c("Chrom","Start","Stop","logSMR","probeID","b_SMR","p_SMR","p_HEIDI")]
# 指定以第4列数据图形中散点纵坐标
data.col <- 4;
# 指定图形在第1个环形生成
track.num <- 1;
# 指定图形在内侧环形生成
side <- "out";
# 指定数据大于4的点以红色显示，小于数据1的点以蓝色显示
by.fold <- 4;
# 绘图
RCircos.Scatter.Plot(smr, data.col,track.num, side, by.fold);

RCircos.Heatmap.Plot(smr, data.col, track.num, side);

RCircos.Histogram.Plot(smr, data.col, track.num, side);





