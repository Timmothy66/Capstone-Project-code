setwd("/Users/rusticbee/Desktop/capstone-Chew/reanaysis/validation with GEO/Brown module/Cluster 4/GSE27887")

library(GEOquery)
my_id <- "GSE27887"
gse <- getGEO(my_id,AnnotGPL = TRUE)

## check how many platforms used
length(gse)
#backupgse<-gse
gse <- gse[[1]]

mypheno_gse <- as.data.frame(pData(gse))
mygene_gse <- as.data.frame(exprs(gse))
fvarLabels(gse) <- make.names(fvarLabels(gse))

# filter samples
mypheno_gse$`treatment:ch1`[which(mypheno_gse$`time:ch1` == "pre-treatment" & mypheno_gse$`tissue:ch1` == "non-lesional skin")]<-"anlpre"
mypheno_gse$`treatment:ch1`[which(mypheno_gse$`time:ch1` == "pre-treatment" & mypheno_gse$`tissue:ch1` == "lesional skin")]<-"alpre"
mypheno_gse$`treatment:ch1`[which(mypheno_gse$`time:ch1` == "post-treatment" & mypheno_gse$`tissue:ch1` == "non-lesional skin")]<-"anlpost"
mypheno_gse$`treatment:ch1`[which(mypheno_gse$`time:ch1` == "post-treatment" & mypheno_gse$`tissue:ch1` == "lesional skin")]<-"alpost"

anlpre<-which(mypheno_gse$`treatment` == "anlpre")
alpre<-which(mypheno_gse$`treatment` == "alpre")
anlpost<-which(mypheno_gse$`treatment` == "anlpost")
alpost<-which(mypheno_gse$`treatment` == "alpost")

gse <- gse[ ,c(anlpre,alpre,anlpost,anlpost)]

qx <- as.numeric(quantile(mygene_gse, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { mygene_gse[which(mygene_gse <= 0)] <- NaN
mygene_gse <- log2(mygene_gse) }

summary(mygene_gse) #log2 transformed
boxplot(mygene_gse,outline=FALSE)

#preparing your data for analysis
gse$group<-factor(mypheno_gse$`treatment:ch1`[c(anlpre,alpre,anlpost,alpost)])
levels(gse$group)<- make.names(c("anlpre","alpre","anlpost","alpost"))

design <- model.matrix(~group + 0, gse)
colnames(design) <- levels(gse$group)

library(limma)
# fit model
fit <- lmFit(gse, design)  

# get top table
cont.matrix <- makeContrasts(contrasts="anlpre-alpre-anlpost-alpost", levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)#empirical bayes
res <- topTable(fit2, adjust="fdr", number=54675)

res <- subset(res, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(res, "gse27887_results.txt", row.names=F, sep="\t")

#整理表达量与形状表
#preparation
hub<-read.csv("/Users/rusticbee/Desktop/capstone-Chew/reanaysis/PPI_STRING/Hub gene/658_AD-atopy/Brown/cluster 4/default_node.csv",header = T)
hub<-read.csv("/Users/rusticbee/Desktop/capstone-Chew/reanaysis/PPI_STRING/functional enrichment/658_AD-atopy/green module/cluster 2/string_node_degrees.tsv",
              header = T,sep = "\t")

hubgeo<-hub$name

reshub<-res[which(res$Gene.symbol %in% hubgeo),]
reshub.fin<-reshub[!duplicated(reshub$Gene.symbol),]
hub$name[-which(hub$name%in%reshub.fin$Gene.symbol)]#找有哪些gene没有
reshub.fin<-as.matrix(reshub.fin)
reshub.fin<-rbind(reshub.fin,res[which(rownames(res)%in%c("207046_at")),])
write.table(reshub.fin, "hub_gse5667_results.txt", row.names=F, sep="\t")

reshub.fin<-read.table("/Users/rusticbee/Desktop/capstone-Chew/reanaysis/validation with GEO/Green module/Cluster 1/GSE22611/hub_gse22611_results.txt",
                       header = T,row.names = 1)
reshub.fin<-as.data.frame(reshub.fin)
#final
finaldat<-t(mygene_gse[which(rownames(mygene_gse) %in% rownames(reshub.fin)),])
finaldat<-as.data.frame(finaldat)
for (i in 1:length(rownames(reshub.fin))) {
  for (j in 1:length(colnames(finaldat))) {
    if(rownames(reshub.fin)[i]==colnames(finaldat)[j]){
      colnames(finaldat)[j]<-reshub.fin$Gene.symbol[i]
             }
  }
}
colnames(finaldat)
colnames(finaldat)[2]<-c("HIST4H4")
finaldat$HIST1H4C<-finaldat$HIST4H4
finaldat$status<-mypheno_gse$`treatment:ch1`



#boxplot
library(ggplot2)
library(ggpubr)
library(magrittr)
library(ggsignif)
library(patchwork)
#install.packages("Rmisc")
library(Rmisc)

finaldat$status<-as.factor(finaldat$status)

p<-list()
for (i in colnames(finaldat)[ncol(finaldat)]){
     for (j in colnames(finaldat)[7:12]) {
       
         
       pl<-list(ggboxplot(finaldat, #数据对象
                                               x = i, # 选择x轴用那一列数据
                                               y = j, #选择y轴用什么数据
                                               fill = 'status', #颜色根据哪一列决定
                                               bxp.errorbar = T, #是否添加error bar
                                               bxp.errorbar.width = 0.2, #error bar的长度
                                               palette = 'npg', #颜色风格
                                               add = 'point' 
                                               )+
        labs(x = "status", # x轴的名字
             y = j # y轴的名字
             )+
        geom_signif(comparisons = list(c('anlpre', 'anlpost'),c('alpre', 'alpost')), # 设置要对比的组
                     step_increase = 0.1,
                     tip_length = c(0.01), #设置显著性那条横线两头向下的长度
                     map_signif_level = F, #设置是否标记显著性的*号，还是直接标记数值
                     test = t.test #设置显著性计算方式
             )+
        theme(
                      plot.title    = element_text(color = 'black', size = 16, hjust = 0.5),
                      plot.subtitle = element_text(color = 'black', size = 13,hjust = 0.5),
                      axis.text.x   = element_text(color = 'black', size = 13, angle = 0),
                      axis.text.y   = element_text(color = 'black', size = 13, angle = 0),
                      axis.title.x  = element_text(color = 'black', size = 13, angle = 0),
                      axis.title.y  = element_text(color = 'black', size = 13, angle = 90),
                      legend.title  = element_text(color = 'black', size  = 15),
                      legend.text   = element_text(color = 'black', size   = 13),
                      axis.line.y = element_line(color = 'black', linetype = 'solid'), # y轴线特征
                      axis.line.x = element_line (color = 'black',linetype = 'solid'), # x轴线特征
                      
              ))
       p<-cbind(p,pl)}
       
 }
Rmisc::multiplot(plotlist = p, layout = matrix(1:6, nrow = 2, byrow = T))





#write.csv(mygene_gse,"megene.csv")












