#setwd
setwd("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/reanaysis/658")
#install WGCNA
#install.packages("WGCNA")
library(WGCNA)

#Get phenotypes of genes
coldata<-read.table("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/AD/WGCNA/for_fanyu/rseq_data/rseq_658samples_geneFPKM_normalized.txt",header = TRUE,sep = "\t")

# log2 transformation
fpkm_matrix <- read.table("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/AD/WGCNA/for_fanyu/rseq_data/rseq_658samples_geneFPKM_normalized.txt",header = T,sep = "\t",row.names = 1)
fpkm_matrix<-as.matrix(fpkm_matrix)
datExpr = t(log2(fpkm_matrix+1))

#coldata<-coldata[which(coldata$AD_status %in% c("case","control")),]
#samp_ind<-which(rownames(datExpr) %in% coldata$RSEQID)
#datExpr<-datExpr[samp_ind,]

gene_mean<-as.vector(apply(datExpr,2,mean))
hist(datExpr,breaks = seq(0,20,0.1))
hist(gene_mean,breaks = seq(0,4140,0.1), xlim = range(0.5))
dev.off()

#remove unexpressed genes. necessary for WGCNA
rm_ind<-which(gene_mean< 0.1 )
datExpr<-datExpr[,-rm_ind]

gene_mad<-as.vector(apply(datExpr,2,mad))#找绝对中位差，用来估计标准差，R返回标准差
rm_ind2<-which(gene_mad==0)
datExpr<-datExpr[,-rm_ind2]

 #Calculate sample distance and cluster the samples
sampleTree=hclust(dist(datExpr),method = "average")
#Plot sample tree
par(cex=0.5)
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 70, col = "red") 
clust= cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
dev.off()


# there are outliers 

# log2 transformation
fpkm_matrix <- read.table("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/AD/WGCNA/for_fanyu/rseq_data/rseq_658samples_geneFPKM_normalized.txt",header = T,sep = "\t",row.names = 1)
fpkm_matrix<-as.matrix(fpkm_matrix)
datExpr = t(log2(fpkm_matrix+1))

#delet outliers, different from other samples
#coldata<-coldata[-c(which(coldata$RSEQID %in% c("MEJ","ME12")))]

#coldata<-coldata[which(coldata$AD_status %in% c("case","atopic control","control")),]
samp_ind<-which(rownames(datExpr) %in% c("MEJ","ME12","ME62","MG74"))
datExpr<-datExpr[-samp_ind,]

gene_mean<-as.vector(apply(datExpr,2,mean))
hist(datExpr,breaks = seq(0,20,0.1))
hist(gene_mean,breaks = seq(0,4140,0.1), xlim = range(0.5))
dev.off()

#remove unexpressed genes. necessary for WGCNA
rm_ind<-which(gene_mean< 0.1 )
datExpr<-datExpr[,-rm_ind]

gene_mad<-as.vector(apply(datExpr,2,mad))
rm_ind2<-which(gene_mad==0)
datExpr<-datExpr[,-rm_ind2]

#choosing soft threshold
powers = c(1:25)#列一堆power为计算相关参数
#列出多SFT power 做无标度拓扑分析，为架构网购找到合适的power
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, corFnc = "bicor",
                        networkType="signed", corOptions=list(use="p", maxPOutliers=0.05))
#列出一个表格，通过mean和Rsq来选择用于转化为无标度网络的power参数
write.table(sft$fitIndices,"658_sft.txt",col.names = T, row.names = F, quote = F, sep = "\t")
par(mfrow = c(1,2))
cex1 = 0.5

#Soft thresholding Plot
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2", type="n", col="red", 
     main=paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col="red" )
abline(h=0.8, col="red")

#Mean connectivity plot
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",
     ylab="Mean Connecivity", type="n", 
     main=paste("Mean Connecivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5], labels = powers, cex = cex1, col="red" )
dev.off()

#Obtain TOM matrix
power=19
#计算拓扑重叠矩阵
TOM = TOMsimilarityFromExpr(datExpr, power = power, verbose = 5, corType = "bicor", 
                            maxPOutliers = 0.05, networkType = "signed")
dissTOM = 1-TOM

# Plot gene Tree
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main="Gene clustering on TOM-based dissimilarity",
     labels=FALSE, hang= 0.04)
dev.off()

#plot gene tree with labels

pamT_mergeF<-cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, 
                           pamRespectsDendro = TRUE, minClusterSize = 30)
#对分层聚类树图的枝修剪作访问
pamT_mergeF_Colors = labels2colors(pamT_mergeF)
#给出颜色
pamT_mergeT0.1<-mergeCloseModules(datExpr, pamT_mergeF_Colors, cutHeight = 0.1)
#合并基因表达网络中特征基因相关性接近的模块（合并）
pamT_mergeT0.1_Colors = pamT_mergeT0.1$colors
pamT_mergeT0.1_MEs = pamT_mergeT0.1$newMEs
#合并之后的信息：颜色，特征基因

plotDendroAndColors(geneTree, cbind(pamT_mergeF_Colors, pamT_mergeT0.1_Colors),
                    c("No Merge","Merge"), dendroLabels = FALSE,
                    hang= 0.3, addGuide = TRUE, guideHang = 0.05)


#Get phenotypes of genes
coldata<-read.csv("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/AR/coldata.csv",header = TRUE)
coldata<-read.table("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/AD/WGCNA/for_fanyu/rseq_data/coldata.txt",header = TRUE,sep = "\t")

coldata$Acne_status<-coldata$Acne
coldata$Acne_status[which(coldata$Acne == "0")]<-"Control"
coldata$Acne_status[which(coldata$Acne == "1")]<-"Case"
#Find correlation of disease
coldata<-coldata[which(coldata$AR_status %in% c("Case", "Control","Intermediate")),]
coldata$AR<-coldata$AR_status
coldata$AR[which(coldata$AR_status %in% c( "Intermediate","Control"))]<-0
coldata$AR[which(coldata$AR_status=="Case")]<-1
colnames(coldata)[which(colnames(coldata)=="RSEQID")]<-"run"

#Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, pamT_mergeT0.1_Colors)$eigengenes

#计算给定数据集中模块的模块特征基因（第一主成分）

MEs=MEs0
dimenME <- dim(MEs)[2]
MEs = orderMEs(MEs0)#重新排序

finaldf <- as.data.frame(matrix(rep("NA",dimenME*2),ncol = 2))
colnames(finaldf) <- c("Module","P-value")
finaldf$Module <- colnames(MEs)

# Calculate pearson correlation coefficients between module eigen-genes and traits
MEs$run<-rownames(MEs)
combined<-merge(MEs,coldata,by="run")
#combined1<-MEs[which(MEs$run %in% coldata$run),]

#combined1<-combined[match(rownames(MEs), combined$run), ]
write.csv(combined,"Acne_658non_MEs.csv",row.names = T)
#moduleTraitCor = WGCNA::cor(MEs,combined$AD, use = "p")#计算协方差，pearson系数,相关性
#moduleTraitPvalue = as.data.frame(corPvalueStudent(moduleTraitCor, nSamples))#对相关性计算Student asymptotic p-value
#moduleTraitCor = as.data.frame(moduleTraitCor)



for(i in 1:dimenME){
  yvar<-combined[,i+1]
  xvar<-combined$AR
  
  mytest <- t.test(yvar~xvar)
  
  
  finaldf[i,2] <- mytest$p.value[1]
  
  
  
  print(i)
}



removeMEchar<-function(x){return(strsplit(x,"ME")[[1]][2])}#拆分
finaldf$Module<-as.vector(sapply(finaldf$Module,removeMEchar))

write.table(finaldf,"Acne_658non_ttest.txt",col.names=T,row.names=F,quote=F,sep="\t")

library(stringr)
splitfun<-function(x){return(str_split(x,"_")[[1]][length(str_split(x,"_")[[1]])])}

mod_name<-finaldf$Module[28]
mod_names<-finaldf$Module[which(finaldf$`P-value`<0.05)]
for(i in 1:length(mod_names)){
  mod_name<-mod_names[i]
  mod_genes<-colnames(datExpr)[which(pamT_mergeT0.1_Colors==mod_name)]
  tmp<-as.vector(sapply(mod_genes,splitfun))
  tmp<-as.data.frame(tmp)
  #write.table(tmp,paste0(mod_name,"pamTmergT_atop_genelist.txt"),row.names = F,col.names = F)
  write.csv(tmp,paste0(mod_name,"pamTmergT_com_genelist.csv"),row.names = F)
}

#FDR
my_pvalue<-read.table("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/reanaysis/reanalysis_wgcna/358/AD_atopy/AD_358atop_ttest.txt")
colnames(my_pvalue)<-my_pvalue[1,]
my_pvalue<-my_pvalue[-1,]

my_fdr<-p.adjust(my_pvalue[,2],method = "fdr")
my_pvalue<-as.data.frame(my_pvalue)
my_pvalue$`P-value`<-my_fdr
colnames(my_pvalue)<-c("Modules","FDR adjusted Pvalue")
write.table(my_pvalue,"AD_358atop_fdr.txt",row.names = F,col.names =T,quote=F,sep="\t")

mod_names<-my_pvalue$Modules[which(my_pvalue$`FDR adjusted Pvalue`<0.05)]

#MEs plot
MEs<-read.csv("/Users/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/reanaysis/reanalysis_wgcna/658/AR_combine/AR_658combine_MEs.csv")
ttest<-read.table("/Usioers/rusticbee/Desktop/bioinfo&biosta/capstone-Chew/reanaysis/reanalysis_wgcna/658/AR_combine/AR_658combine_fdr.txt",sep = "\t")

removeMEchar<-function(x){return(strsplit(x,"ME")[[1]][2])}#拆分
colnames(MEs)<-as.vector(sapply(colnames(MEs),removeMEchar))

ttest<-ttest[-1,]
mod_names<-ttest$V1[which(ttest$V2<0.05)]
MEs_sign<-MEs[,which(colnames(MEs) %in% mod_names)]
rownames(MEs_sign)<-MEs[,2]
MEs_sign$AR<-MEs[,33]
MEs_sign$AR_sta<-MEs[,32]
MEs_sign$AR_sta[which(MEs_sign$AR==0)]<-"Combine"
MEs_sign<-MEs_sign[order(MEs_sign$AR),]

ggplot(MEs,aes(x=AD_sta,y=MEskyblue3,fill=AD_sta))+theme_bw()+
  geom_violin()+
  geom_signif(comparisons=list(c("atopic control","case")),
              step_increase = 0.1,
              test="t.test",
              map_signif_level=F)

for (i in colnames(MEs_sign)[-c(16,17)]) {
  p<-ggplot(MEs_sign,aes(x=AR_sta,y=eval(parse(text = i)),fill=AR_sta))+theme_bw()+
    ylab(i)+
    geom_violin()+
    geom_signif(comparisons=list(c("Combine","Case")),
                step_increase = 0.1,
                test="t.test",
                map_signif_level=F)
  ggsave(paste0("./",i,'.png'),p,width = 8,height = 7)
}
