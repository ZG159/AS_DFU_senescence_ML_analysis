
library(Seurat)
## 01绘图SASP基因表达------------------------------


VlnPlot(sce.all_int, features = target_gene_DFU_IPH,ncol = 2) #这里我们做DFU组的图
VlnPlot(sce.all_int_Con, features = target_gene_DFU_IPH,ncol = 2) #这里我们做Con组的图


library(ggplot2)
# 提取表达矩阵
data001 <- FetchData(sce.all_int, vars = c(target_gene_DFU_IPH, "celltype", "group"))

data001_1 <- FetchData(sce.all_int, vars = c("CALCOCO1", "celltype", "group"))
data001_2 <- FetchData(sce.all_int, vars = c("CCND2", "celltype", "group"))
data001_3 <- FetchData(sce.all_int, vars = c("FUCA1", "celltype", "group"))




#绘制 CALCOCO1-------------------------

ggplot(data001_1, aes(x = celltype, y = CALCOCO1, fill = group)) +
  geom_boxplot(
    position = position_dodge(0.8),  # 分组间距
    width = 0.7,                     # 箱线图宽度
    alpha = 0.8,                      # 透明度
    outlier.shape = NA                # 隐藏异常值点
  ) +
  
  stat_compare_means(
    aes(group = group),               # 指定分组变量
    method = "wilcox.test",           # 使用Wilcoxon检验
    label = "p.signif",               # 显示显著性符号(*)
    label.y = max(meta_data_sensore$score) * 2.5,  # 标签位置
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 自定义颜色
 
  coord_cartesian(ylim = c(0, 0.75)) + # 改变纵坐标轴范围
  labs(
    x = "celltype",
    y = "CALCOCO1 Expression value",
    fill = "group"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

target_gene_DFU_IPH


#绘制 CCND2-------------------------

ggplot(data001_2, aes(x = celltype, y = CCND2, fill = group)) +
  geom_boxplot(
    position = position_dodge(0.8),  # 分组间距
    width = 0.7,                     # 箱线图宽度
    alpha = 0.8,                      # 透明度
    outlier.shape = NA                # 隐藏异常值点
  ) +
  
  stat_compare_means(
    aes(group = group),               # 指定分组变量
    method = "wilcox.test",           # 使用Wilcoxon检验
    label = "p.signif",               # 显示显著性符号(*)
    label.y = max(meta_data_sensore$score) * 6.8,  # 标签位置
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 自定义颜色
  
  coord_cartesian(ylim = c(0, 2)) + # 改变纵坐标轴范围
  labs(
    x = "celltype",
    y = "CCND2 Expression value",
    fill = "group"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )





#绘制 FUCA1-------------------------

ggplot(data001_3, aes(x = celltype, y = FUCA1, fill = group)) +
  geom_boxplot(
    position = position_dodge(0.8),  # 分组间距
    width = 0.7,                     # 箱线图宽度
    alpha = 0.8,                      # 透明度
    outlier.shape = NA                # 隐藏异常值点
  ) +
  
  stat_compare_means(
    aes(group = group),               # 指定分组变量
    method = "wilcox.test",           # 使用Wilcoxon检验
    label = "p.signif",               # 显示显著性符号(*)
    label.y = max(meta_data_sensore$score) * 6.8,  # 标签位置
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +  # 自定义颜色
  
  coord_cartesian(ylim = c(0, 2)) + # 改变纵坐标轴范围
  labs(
    x = "celltype",
    y = "FUCA1 Expression value",
    fill = "group"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  ) 


# 绘制不分面小提琴图（不区分细胞类型）------------------------

library(ggpubr)
ggplot(data_long, aes(x = gene, y = value, fill = group)) +
  
  geom_violin(scale = "width", trim = TRUE) +
  scale_fill_manual(values = c("#1F77B4", "#FF7F0E")) +
  theme_minimal() +
  stat_compare_means(
    aes(group = group), 
    method = "wilcox.test", 
    label = "p.signif",  # 显示星号
    label.y = 4        # 位置调整
  )+
  coord_cartesian(ylim = c(0, 4))
