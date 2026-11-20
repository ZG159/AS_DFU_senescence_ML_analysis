
library(ggpubr)
library(dplyr)
#1 计算exp_134的差异表达图-------------------------------

data_make_comparation <- expr_134[target_gene,] %>% 
  t() %>% 
  as.data.frame()


data_make_comparation$group <-  c(rep("Con",10),rep("DFU",11))


plot_target_gene <- function(target_gene){
  
  p <- ggboxplot(data_make_comparation, x = "group", y = target_gene, width = 0.6, 
                 color = "black",#轮廓颜色
                 fill="group",#填充
                 palette =c("#E7B800", "#00AFBB"),#分组着色
                 xlab = F, #不显示x轴的标签
                 bxp.errorbar=T,#显示误差条
                 bxp.errorbar.width=0.5, #误差条大小
                 size=1, #箱型图边线的粗细
                 outlier.shape=NA, #不显示outlier
                 legend = "right"
  ) #图例放右边
  
  my_comparisons <- list(c("DFU","Con"))
  
  p
  p + stat_compare_means(comparisons = my_comparisons,label = "p.signif")
  
}
target_gene
plot_target_gene("CCND2")
plot_target_gene("FUCA1")
plot_target_gene("CALCOCO1")

