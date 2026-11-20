library(ggplot2)
library(pROC)


dd_163 <- target_gene_data1
dd <- dd_163

dd$group <- as.factor(dd$group)

dd_fit<- glm(group ~ CALCOCO1 + FUCA1 + CCND2 ,data = dd,family = binomial)
dd$Prediction_model <- predict(dd_fit, type="response")


roc_obj <- roc(dd$group, dd$Prediction_model)




as.numeric(auc(roc_obj))


roc_obj_smooth <- smooth(roc_obj)

roc_obj_smooth
auc_value2 <- as.numeric(auc(roc_obj_smooth2))
# 提取 ROC 曲线的坐标
roc_data <- data.frame(
  specificity = 1 - roc_obj_smooth$specificities,
  sensitivity = roc_obj_smooth$sensitivities
)

roc_data 
ggplot(roc_data, aes(x = specificity, y = sensitivity)) +
  geom_line(color = "#4e91ff", linewidth = 1) +  # ROC 曲线
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 对角线（随机分类器）
  labs(
    title = "predict model",
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
