setwd("/Users/rusticbee/Desktop/capstone-Chew/reanaysis/Comparison/AD_Atopy/Black module/Cluster 3/non")
#导入原始数据
fpkm<-read.table("/Users/rusticbee/Desktop/capstone-Chew/AD/WGCNA/for_fanyu/rseq_data/rseq_658samples_geneFPKM_normalized.txt",header = T,sep="\t",row.names = 1)
coldata<-read.table("/Users/rusticbee/Desktop/capstone-Chew/AD/WGCNA/for_fanyu/rseq_data/coldata.txt",header=T,sep="\t")
hub<-read.csv("/Users/rusticbee/Desktop/capstone-Chew/reanaysis/PPI_STRING/Hub gene/658_AD-atopy/Black/cluster 3/default_node.csv",header = T)

#整理数据
#整理fpkm
library(stringr)
splitfun<-function(x){
  return(str_split(x,"_")[[1]][length(str_split(x,"_")[[1]])])
}
fpkm.matrix<-t(as.matrix(fpkm))
colnames(fpkm.matrix)<-as.vector(sapply(colnames(fpkm.matrix),splitfun))

hubge<-hub$name
fpkm.chos<-fpkm.matrix[,which(colnames(fpkm.matrix) %in% hubge)]
genemean<-as.vector(apply(fpkm.chos,2,mean))
fpkm.chose<-as.data.frame(fpkm.chos[,-which(genemean==0)])
fpkm.chose<-as.data.frame(fpkm.chos)

coldata<-coldata[which(coldata$AD_status %in% c("case","control")),]
fpkm.chose1<-fpkm.chose[which(rownames(fpkm.chose) %in% coldata$RSEQID),]
fpkm.chose1$AD_status<-coldata$AD_status

#boxplot
#library(ggplot2)
#library(ggpubr)
#library(magrittr)
#library(ggsignif)

#install.packages("pacman")
#pacman::p_load(ggplot2, ggpubr, magrittr, ggsignif)

fpkm.chose1$AD_status<-as.factor(fpkm.chose1$AD_status)

for (i in colnames(fpkm.chose1)[-ncol(fpkm.chose1)]){
  p<-ggboxplot(fpkm.chose1, #数据对象
               x = 'AD_status', # 选择x轴用那一列数据
               y = 'fpkm.chose1[[i]]', #选择y轴用什么数据
               fill = 'AD_status', #颜色根据哪一列决定
               bxp.errorbar = T, #是否添加error bar
               bxp.errorbar.width = 0.2, #error bar的长度
               palette = 'npg', #颜色风格
               add = 'point' # 是否添加boxplot上面的点点
  )+
    labs(title = 'Gene expression in two conditions', # 添加主标题
         subtitle = 'Plot of GE by condition', # 添加次标记
         x = 'AD status', # x轴的名字
         y = 'Gene expression' # y轴的名字
    )+
    geom_signif(comparisons = list(c('case', 'control')), # 设置要对比的组
                #y_position = c(400,402,404), #设置3个显著性标记的高度
                step_increase = 0.1,
                tip_length = c(0.01), #设置显著性那条横线两头向下的长度
                map_signif_level = F, #设置是否标记显著性的*号，还是直接标记数值
                test = t.test #设置显著性计算方式
    )+
    theme(
      plot.title    = element_text(color = 'black', size   = 16, hjust = 0.5),
      plot.subtitle = element_text(color = 'black', size   = 13,hjust = 0.5),
      axis.text.x   = element_text(color = 'black', size = 13, angle = 0),
      axis.text.y   = element_text(color = 'black', size = 13, angle = 0),
      axis.title.x  = element_text(color = 'black', size = 13, angle = 0),
      axis.title.y  = element_text(color = 'black', size = 13, angle = 90),
      legend.title  = element_text(color = 'black', size  = 15),
      legend.text   = element_text(color = 'black', size   = 13),
      axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
      axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
      panel.border = element_rect(linetype = 'solid', size = 1.2,fill = NA) # 图四周框起来
    )
  ggsave(paste0('./',i,'.png'),p,width = 8,height = 7)
}
