library(COSG)
library(harmony)
library(Seurat)
library(dplyr)
library(stringr)
library(patchwork)
library(tidyr)
library(ggsci)
library(cowplot)
library(data.table)
library(ggplot2)
library(future)
library(harmony)
library(cowplot)
library(RhpcBLASctl)


#本文参考“生信漫漫学”的“单细胞亚群注释三种方法汇总”（https://mp.weixin.qq.com/s/1JRF5aYAkR0U8z2I1nawgQ）


# 02 在线网站注释--------------------------------------------

## a.利用COSG函数获取每一簇的TOP基因-------------------------------------------------
sce.all_int[["RNA"]] <- JoinLayers(sce.all_int[["RNA"]])#这一歩不进行，下面差异分析会报错

marker_cosg <- cosg(
  sce.all_int,
  groups='all',## 默认是在active.idents中的全部，也可以设置多组比较 e.g. c('A cells', 'B cells', 'C cells')
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100)

## Top6 genes
top_6 <- unique(as.character(apply(marker_cosg$names,2,head,6)))
# width <-0.006*dim(sce.all.int)[2];width
# height <- 0.25*length(top_10)+4.5;height

#是否需要画热图？有待进一步整理？
# width <- 15+0.5*length(unique(Idents(sce.all.int)));width
# height <- 8+0.1*length(top_6);height
# sce.Scale <- ScaleData(sce.all.int ,features =  top_6  ) #为了画热图，对数据进行标准化 
# DoHeatmap(  sce.Scale , top_6 , 
#             size=3)

# ggsave(filename=paste0(pro,'DoHeatmap_check_top6_markers_by_clusters.pdf') ,
#        # limitsize = FALSE,
#        units = "cm",width=width,height=height)

## b.根据top基因在ACT网站在线注释-----------------------------------------------

library(dplyr)
cat(paste0('cluster',0:20,':',
           unlist(apply(marker_cosg$names,2,function(x){
             paste(head(x),collapse=',')
           })),'\n'))
#以上进入ACT网站，进行注释得到结果文件"任务：DFU and Cell Senescence/resource/ACT 注释结果.txt"
#下面我们来读取文件：

ACT_annotation <- read_excel("resource/ACT注释结果（已删除NA）.xlsx")


