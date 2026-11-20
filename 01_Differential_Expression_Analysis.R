
library(GEOquery)
library(readxl)  
library(dplyr)
library(tibble)
##GSE163154数据下载及差异基因分析----------------------
expr_163<-getGEO("GSE163154",destdir = '.',AnnotGPL = T,getGPL = T)
pdata1 <- pData(expr_163[[1]])#提取临床信息，包含分组信息
fdata1<-fData(expr_163[[1]])#提取平台信息,其中包含芯片探针的注释信息
expr_163<-exprs(expr_163[[1]])%>%as.data.frame()#提取表达矩阵
expr_163$ID <- rownames(expr_163)
fdata2 <- fdata1[,c(1,3)]

expr_163_1 <- expr_163 %>% 
  inner_join(fdata2,by="ID")
expr_163_1 <- expr_163_1 %>% 
  select(44:45,1:43)

rownames(expr_163_1) <- expr_163_1$ID

expr_163_1_1 <-  expr_163_1 %>% 
  mutate(rowMean =rowMeans(.[,-c(1:2)])) %>% 
  ## 把表达量的平均值按从大到小排序
  arrange(desc(rowMean)) %>% 
  ## 去重，symbol留下第一个
  distinct(`Gene symbol`,.keep_all = T) %>% 
  ## 反向选择去除rowMean这一列
  select(-rowMean) 

expr_163_1_1 <- expr_163_1_1 %>% 
  filter(ID != "ILMN_1683271")



expr_163_1_1 <- expr_163_1_1 %>% 
  select(- ID) 

expr_163_1_1 <- as.data.frame(expr_163_1_1)
rownames(expr_163_1_1) <- expr_163_1_1$`Gene symbol`
expr_163_1_1 <- expr_163_1_1[,-1]
saveRDS(expr_163_1_1,file = "../任务：DFU+IPH/expr_163.Rdata")
## 自动log化，解释
## 注意此时数据的格式是矩阵
ex1 <- expr_163_1_1
ex1 <- as.matrix(ex1)
qx1 <- as.numeric(quantile(ex1, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx1[5] > 100) ||
  (qx1[6]-qx1[1] > 50 && qx1[2] > 0) ||
  (qx1[2] > 0 && qx1[2] < 1 && qx1[4] > 1 && qx1[4] < 2)

## 开始判断
if (LogC) { 
  ex1[which(ex1 <= 0)] <- NaN
  ## 取log2
  exprSet1 <- log2(ex1)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")  
}

library(limma) 
boxplot(ex1,outline=FALSE, notch=T, las=2)
### 该函数默认使用quntile 矫正差异 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
## 这步把矩阵转换为数据框很重要
exprSet1 <- as.data.frame(ex1)


group1 <- c(rep("Con",16),rep("IPH",27))
group1
group1 <- factor(group1,levels = c("Con","IPH")) 
group1

design1<-model.matrix(~0+group1)
colnames(design1)<-levels(group1)
design1
contrast.matrix1<-makeContrasts(IPH-Con,levels=design1)
contrast.matrix1

fit1<-lmFit(exprSet1,design1)#线性模型拟合
fit1<-contrasts.fit(fit1,contrast.matrix1)
fit1<-eBayes(fit1)#贝叶斯检验

#最终得到差异分析结果
allDiff1=topTable(fit1,adjust='fdr',coef=1,number=Inf)
diffgene_diff_all1 <-  allDiff1 %>% 
  filter(AveExpr>0)

diffgene_diff_all1$group <- ifelse(diffgene_diff_all1$P.Value <0.05 & abs(diffgene_diff_all1$logFC) > 0.5,    
                                   ifelse(diffgene_diff_all$logFC>1,"Up","Down"),
                                   "Not")


table(diffgene_diff_all1$group)

diffgene1_up <- diffgene_diff_all1[diffgene_diff_all1$group=="Up",]
diffgene1_up <- na.omit(diffgene1_up)
diffgene_up1 <- rownames(diffgene1_up)


diffgene_down1 <- diffgene_diff_all1[diffgene_diff_all1$group=="Down",]
diffgene_down1 <- na.omit(diffgene_down1)
diffgene_down1 <- rownames(diffgene_down1)

#得到差异基因
diffgene1 <- c(diffgene_up1,diffgene_down1)


## GSE134431数据下载及差异基因分析------------------------------
expr_134 <- read_excel("resource/GSE134431_181016_diabetic_rna-seq_results.gene.rpkm.xlsx") 

expr_134 <- expr_134[-1,]

library(janitor)
expr_134 <- row_to_names(expr_134, row_number = 1)  # 将第一行作为列名
expr_134_type <- expr_134 %>% 
  select(1,23,24,25)

expr_134 <- expr_134 %>% 
  select(-c(23,24,25))

expr_134 <- expr_134 %>% 
  tibble::column_to_rownames("ID")
expr_134 <- expr_134 %>% 
  dplyr::mutate(across(everything(),as.numeric))


##差异基因分析
metadata2 <- data.frame(sample= colnames(expr_134),group=c(rep("DFU",7),rep("Control",13),"DFU"))
metadata2$group <- factor(metadata2$group,levels = c("Control","DFU")) 

## 1.对dds进行处理，并看数据分布
dds2 <-DESeqDataSetFromMatrix(countData=round(expr_134), 
                              colData=metadata2, 
                              design=~group, #注意design的对象是什么，不是默认的
                              tidy=FALSE)#当第一列为基因时，tidy=true



dds2 <- dds2[rowSums(counts(dds2))>1,] #筛选符合条件的数据
nrow(dds2)

## result得到差异结果
dds2 <- DESeq(dds2)
contrast <- c("group","DFU","Control")#依次是，1.分组信息(metadata中的列) 2.处理组，3.对照组
dd2 <- DESeq2::results(dds2, contrast=contrast, alpha = 0.05)
dd2<- as.data.frame(dd2)
head(dd2)

## 设置阈值，定义差异基因
logFC_cutoff <- 2
dd2$change <-  as.factor(ifelse(dd2$pvalue < 0.05 & abs(dd2$log2FoldChange) > logFC_cutoff,
                                ifelse(dd2$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
table(dd2$change)
dd2_2 <- dd2
dd2_2$change <-  as.factor(ifelse(dd2_2$pvalue < 0.05 & abs(dd2_2$log2FoldChange) > 0.5,
                                  ifelse(dd2_2$log2FoldChange > 0.5 ,'UP','DOWN'),'NOT'))
table(dd2_2$change)
#差异基因
diffgene2 <- rownames(dd2_2[dd2_2$change!="NOT",]) 


#得到衰老基因集
SCG1 <- read.csv("Search Signatures of Cellular Senescence.csv")
SCG1 <- SCG1$Symbol
length(SCG1)

#取基因交集
gene <- intersect(SCG1,intersect(diffgene1,diffgene2))

