#调用pROC包；
library(pROC)
#调用ggplot2;
library(ggplot2)



# 数据准备
target_gene_data <- expr_134[target_gene,]
target_gene_data <- target_gene_data %>% 
  t() %>% 
  as.data.frame()
expr_134_group <- c(rep("DFU",7),rep("Control",13),"DFU")
expr_134_group <- factor(expr_134_group,levels = c("Control","DFU")) 
target_gene_data$group <- expr_134_group
#查看数据类型；
class(target_gene_data$group)
#[1] "factor"
target_gene

# 1.CALCOCO1 -------------------------------------------------------------

rocobj1 <- roc(target_gene_data$group,target_gene_data$CALCOCO1,ci=T,auc=T)
rocobj1


as.numeric(auc(rocobj1))

roc_obj_smooth1 <- smooth(rocobj1)

roc_obj_smooth1
auc_value1 <- as.numeric(auc(roc_obj_smooth1))
# 提取 ROC 曲线的坐标
roc_data1 <- data.frame(
  specificity = 1 - roc_obj_smooth1$specificities,
  sensitivity = roc_obj_smooth1$sensitivities
)

roc_data1 
ggplot(roc_data1, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#4e91ff", linewidth = 1) +  # ROC 曲线
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 对角线（随机分类器）
  labs(
    title = "CALCOCO1",
    x = "1 - Specificity ",
    y = "Sensitivity "
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 10), #设置坐标轴的字体大小
    axis.text.y = element_text(size = 10),
    
    axis.title.x = element_text(size = 10),# 修改坐标轴标题的字体大小
    axis.title.y = element_text(size = 10),
    
    plot.title = element_text(hjust = 0.5),
    # 添加边框
    panel.border = element_rect(color = "black", size = 0.5, fill = NA),  # 显式添加边框
    axis.line = element_blank()  # 移除坐标轴线
  ) +
  # 添加 AUC 和 95% 置信区间文本
  annotate(
    "text",
    x = 0.7,  # 文本的 x 位置
    y = 0.3,  # 文本的 y 位置
    label = paste("AUC:", round(auc_value1, 3),
                  "\n95% CI:", 
                  round(ci.auc(rocobj1)[1], 3), "-",
                  round(ci.auc(rocobj1)[3], 3)),
    size = 4,
    color = "black",
    hjust = 0
  )

# 2.FUCA1-------------------------------------------------------------

rocobj2 <- roc(target_gene_data$group,target_gene_data$FUCA1,ci=T,auc=T)
rocobj2
as.numeric(auc(rocobj2))

roc_obj_smooth2 <- smooth(rocobj2)

roc_obj_smooth2
auc_value2 <- as.numeric(auc(roc_obj_smooth2))
# 提取 ROC 曲线的坐标
roc_data2 <- data.frame(
  specificity = 1 - roc_obj_smooth2$specificities,
  sensitivity = roc_obj_smooth2$sensitivities
)

roc_data2 
ggplot(roc_data2, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#4e91ff", linewidth = 1) +  # ROC 曲线
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 对角线（随机分类器）
  labs(
    title = "FUCA1",
    x = "1 - Specificity ",
    y = "Sensitivity "
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 10), #设置坐标轴的字体大小
    axis.text.y = element_text(size = 10),
    
    axis.title.x = element_text(size = 10),# 修改坐标轴标题的字体大小
    axis.title.y = element_text(size = 10),
    
    plot.title = element_text(hjust = 0.5),
    # 添加边框
    panel.border = element_rect(color = "black", size = 0.5, fill = NA),  # 显式添加边框
    axis.line = element_blank()  # 移除坐标轴线
  ) +
  # 添加 AUC 和 95% 置信区间文本
  annotate(
    "text",
    x = 0.7,  # 文本的 x 位置
    y = 0.3,  # 文本的 y 位置
    label = paste("AUC:", round(auc_value2, 3),
                  "\n95% CI:", 
                  round(ci.auc(rocobj2)[1], 3), "-",
                  round(ci.auc(rocobj2)[3], 3)),
    size = 4,
    color = "black",
    hjust = 0
  )

# 2.FUCA1-------------------------------------------------------------

rocobj3 <- roc(target_gene_data$group,target_gene_data$CCND2,ci=T,auc=T)
rocobj3
as.numeric(auc(rocobj3))

roc_obj_smooth3 <- smooth(rocobj3)

roc_obj_smooth3
auc_value3 <- as.numeric(auc(roc_obj_smooth3))
# 提取 ROC 曲线的坐标
roc_data3 <- data.frame(
  specificity = 1 - roc_obj_smooth3$specificities,
  sensitivity = roc_obj_smooth3$sensitivities
)

roc_data3 
ggplot(roc_data3, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#4e91ff", linewidth = 1) +  # ROC 曲线
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 对角线（随机分类器）
  labs(
    title = "CCND2",
    x = "1 - Specificity ",
    y = "Sensitivity "
  ) +       
  theme_minimal() +      
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 10), #设置坐标轴的字体大小
    axis.text.y = element_text(size = 10),
    
    axis.title.x = element_text(size = 10),# 修改坐标轴标题的字体大小
    axis.title.y = element_text(size = 10),
    
    plot.title = element_text(hjust = 0.5),
    # 添加边框
    panel.border = element_rect(color = "black", size = 0.5, fill = NA),  # 显式添加边框
    axis.line = element_blank(),  # 移除坐标轴线
    panel.grid.major = element_blank(), # 新增：去除主网格线
    panel.grid.minor = element_blank(), # 新增：去除次网格线
  ) +
  # 添加 AUC 和 95% 置信区间文本
  annotate(
    "text",
    x = 0.7,  # 文本的 x 位置
    y = 0.3,  # 文本的 y 位置
    label = paste("AUC:", round(auc_value3, 3),
                  "\n95% CI:", 
                  round(ci.auc(rocobj3)[1], 3), "-",
                  round(ci.auc(rocobj3)[3], 3)),
    size = 4,
    color = "black",
    hjust = 0
  )




