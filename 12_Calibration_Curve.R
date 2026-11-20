# 01校准曲线 134-----------------------------------

library(riskRegression)

dd_Calibration <- data_make_comparation
dd_Calibration$group <- as.factor(dd_Calibration$group)


Calibration_fit <- glm(group ~ CALCOCO1 + FUCA1 + CCND2,
                       data = dd_Calibration, family = binomial)

Calibration_fit_score <- Score(list("fit"= Calibration_fit),
                               formula = group ~ 1,
                               data = dd_Calibration,
                               metrics = c("auc","brier"),
                               #summary = c("risks","IPA","riskQuantile","ibs"),
                               plots = "calibration",
                               null.model = T,
                               conf.int = T,
                               B = 500
)



plotdata <- plotCalibration(Calibration_fit_score)

library(ggplot2)

ggplot(plotdata$plotFrames$fit, aes(x=Pred,y=Obs))+
  geom_line(color="tomato",linewidth=1.5)+
  scale_x_continuous(limits = c(0,1),name = "Predicted Risk")+
  scale_y_continuous(limits = c(0,1),name = "Observerd Risk")+
  geom_abline(slope = 1,intercept = 0,lty=2)+
  geom_rug(color="grey")+
  theme_bw()

target_gene

#另一种画法：

library(rms)
# 加载rms包
Model <- lrm(group ~ CALCOCO1 + FUCA1 + CCND2, x = T, y = T, data = dd_Calibration) #构建模型Model 
#输出模型结果
cal<-calibrate(Model,data = dd_Calibration)
plot(cal,  xlim=c(0,1.0),ylim=c(0,1.0),  xlab = "Predicted Probability", ylab = "Observed Probability") #输出校准曲线图片

# 02校准曲线 163-----------------------------------（是以136列线图的模型做的）


dd_Calibration_163 <- dd_163
dd_Calibration_163$group <- as.factor(dd_Calibration_163$group)


library(rms)
# 加载rms包
Model_163 <- lrm(group ~ CALCOCO1 + FUCA1 + CCND2, x = T, y = T, data = dd_Calibration_163) #构建模型Model 
#输出模型结果
cal_163<-calibrate(Model,data = dd_Calibration_163)
plot(cal_163,  xlim=c(0,1.0),ylim=c(0,1.0),  xlab = "Predicted Probability", ylab = "Observed Probability") #输出校准曲线图片



