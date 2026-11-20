library(ggplot2)
library(pROC)


dd_136 <- target_gene_data
dd <- dd_136
dd$group <- as.factor(dd$group)

dd_fit<- glm(group ~ CALCOCO1 + FUCA1 + CCND2 ,data = dd,family = binomial)
dd$Prediction_model <- predict(dd_fit, type="response")


roc_obj <- roc(dd$group, dd$Prediction_model)


roc_data <- data.frame(
  sensitivity = roc_obj$sensitivities,
  specificity = roc_obj$specificities
)
roc_data$fpr <- 1 - roc_data$specificity # 计算FPR

auc_value <- as.numeric(auc(roc_obj))   # 计算AUC
 ci_value <- ci(roc_obj)   #获取置信区间

#使用ggplot2绘制，并使用geom_smooth进行平滑



ggplot(roc_data, aes(x = fpr, y = sensitivity)) +
  # 添加对角线（红色虚线）
 
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # 对角线（随机分类器）
  # 添加平滑的ROC曲线（蓝色）
  geom_smooth(method = "loess", se = FALSE, span = 0.8, 
              color = "#4e91ff", size = 1.2, alpha = 0.8) +
  # 设置坐标轴范围
  xlim(0, 1) +
  ylim(0, 1) +
  # 添加标签和标题
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    title = "predict model"
  ) +
  # 添加AUC和置信区间标注
  annotate("text", x = 0.6, y = 0.2, 
           label = paste0("AUC: ", round(auc_value, 2)), 
           size = 4, hjust = 0) +
  annotate("text", x = 0.6, y = 0.1, 
           label = paste0("95% CI: ", round(ci_value[1], 3), " - ", round(ci_value[3], 3)), 
           size = 4, hjust = 0) +
  # 使用简洁主题并添加网格
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
    plot.background = element_rect(fill = "white", colour = "lightgray", size = 0.7)

)













