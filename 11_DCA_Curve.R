
# 原文链接：https://blog.csdn.net/zfyyzhys/article/details/144006922
library(dcurves)

# DCA图136----------------------------------
dd_dca <- data_make_comparation
dd_dca$group <- as.factor(dd_dca$group)

dd_dca_fit<- glm(group ~ CALCOCO1 + FUCA1 + CCND2 ,data = dd_dca,family = binomial)
dd_dca$Prediction_model <- predict(dd_dca_fit, type="response")

dca(group ~ Prediction_model ,
    data = dd_dca,
    thresholds = seq(0, 1.0, by = 0.01),# 设定不同的阈值
    label = list(model_train = "Train dataset")) %>%
  plot(smooth = TRUE)



# DCA图163----------------------------------

library(dplyr)
dd_dca_163 <- dd_163
dd_dca_163$group <- as.factor(dd_dca_163$group)

dd_dca_fit_163<- glm(group ~ CALCOCO1 + FUCA1 + CCND2 ,data = dd_dca_163,family = binomial)
dd_dca_163$Prediction_model <- predict(dd_dca_fit_163, type="response")

dca(group ~ Prediction_model ,
    data = dd_dca_163,
    thresholds = seq(0, 1.0, by = 0.01),# 设定不同的阈值
    label = list(model_train = "Train dataset")) %>%
  plot(smooth = TRUE)



