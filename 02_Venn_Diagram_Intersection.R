library(ggplot2)
install.packages("ggvenn")
library(ggvenn)
# 准备数据为命名列表
venn_data <- list(
  SVM = svm_geneids,
  Lasso = lasso_geneids,
  RF = rf_geneids
)

# 绘制韦恩图
ggvenn(
  venn_data,
  fill_color = c("#1b9e77", "#d95f02", "#7570b3"), # 自定义颜色
  stroke_size = 0.8,      # 圆圈线宽
  set_name_size = 4.5,    # 集合名称大小
  show_percentage = FALSE  # 不显示百分比
) 

target_gene
