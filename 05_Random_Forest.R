
#完整代码：
#https://mp.weixin.qq.com/s?__biz=MzkyNDc1OTYzMg==&mid=2247489285&idx=1&sn=d4881050cdc6bc6d49e7e9065110170b&chksm=c0ada20eee5e0cdc3a31e20b723cfd551156663742d9f5fd0c397157e4fc3901c9b24d562dc9#rd

set.seed(2025)
devtools::install_github("ModelOriented/randomForestExplainer")
library(randomForest)
library(randomForestExplainer)
library(dplyr)
library(ggplot2)
#数据准备
exp_rf <- exp_lasso
exp_rf$group <- as.factor(exp_lasso$group)
names(exp_rf) <- gsub("-", "_", names(exp_rf))


# 参数优化函数（寻找最佳特征子集大小）
tune_params <- function(exp_svm) {
  set.seed(2025)
  tuneRF(
    x = x_all,  # 特征矩阵
    y = exp_rf$group,   # 分组标签
    ntreeTry = 100,  # 初始树数量
    stepFactor = 1.5, # 每次增加50%的mtry
    improve = 0.05,   # 误差改善<5%时停止
    trace = TRUE
  )
}

tune_params(exp_rf) 
# 获取最优mtry值
best_mtry <- tune_params(exp_rf) %>% 
  as.data.frame() %>% 
  dplyr::slice(which.min(OOBError))%>% 
  pull(mtry) 

# 构建最终模型
rf_model <- randomForest(
  group ~ .,  # 使用所有特征预测分组
  data = exp_rf,
  ntree = 500,        # 建造500棵决策树
  mtry = best_mtry,   # 使用优化后的特征数
  importance = TRUE,  # 计算特征重要性
  do.trace = 50       # 每50棵树打印进度
)

importance(rf_model)
# 设置重要性阈值：中位数 
imp_scores <- importance(rf_model)[, "MeanDecreaseGini"]
threshold <- median(imp_scores) 
rf_geneids <- names(imp_scores[imp_scores > threshold])
rf_geneids

# OOB误差曲线（展示模型收敛过程）
oob_error <- rf_model$err.rate[, "OOB"]
err_df <- data.frame(Trees = 1:100, OOB_Error = oob_error)
p_err <- ggplot(err_df, aes(Trees,  OOB_Error)) +
  geom_line(color = "black", linewidth=0.7) +
  geom_vline(xintercept=which.min(err_df$Error), linetype=2)
p_err

plot(rf_model, main="OOB Error Rate vs. Number of Trees")

# 特征重要性图谱
imp_df <- data.frame(Gene=names(imp_scores), 
                     MeanDecreaseGini=imp_scores)
imp_df_1 <- imp_df %>% 
  filter(rownames(imp_df) %in% rf_geneids)

ggplot(imp_df_1, aes(MeanDecreaseGini, reorder(Gene, MeanDecreaseGini))) +
  geom_segment(aes(xend=0, yend=Gene), color="grey70") +
  geom_point(aes(color=MeanDecreaseGini), size=4, alpha=0.8) +
  scale_color_gradient(low="#00AFBB", high="#FC4E07")

ggplot(imp_df_1, aes(MeanDecreaseGini, reorder(Gene, MeanDecreaseGini))) +
  geom_segment(aes(xend=0, yend=Gene), color="grey70") +
  geom_point(aes(color=MeanDecreaseGini), size=4, alpha=0.8) +
  scale_color_gradient(low="#00AFBB", high="#FC4E07") +
  labs(x = "Mean Decrease Gini", y = "Gene") +  # 添加横纵坐标标签
  theme(
    panel.background = element_rect(fill = "white"),  # 设置背景为纯白色
    panel.grid.major = element_blank(),  # 移除网格线
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),  # 添加坐标轴线
    axis.text = element_text(color = "black"),  # 设置坐标轴文本颜色
    axis.title = element_text(color = "black")  # 设置坐标轴标题颜色
  )


intersect(rf_geneids,lasso_geneids)






