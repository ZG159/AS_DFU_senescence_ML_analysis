#参考资料：https://ayueme.github.io/R_clinical_model/nomogram-logistic.html

# 01绘制列线图136-------------------------------
dd_nomo <- data_make_comparation
dd_nomo$group <- as.factor(dd_nomo$group)

library(rms)

rms_fit <- rms::lrm(group ~ CALCOCO1 + FUCA1 + CCND2 ,
                    data = dd_nomo)
#fit1
dd_nomo <- datadist(dd_nomo)
options(datadist="dd_nomo")
nom <- nomogram(rms_fit, fun=plogis,
                fun.at=c(0.001,0.1,0.25,0.5,0.75,0.9,0.99),
                lp=T, # 是否显示线性预测值
                maxscale = 100, # 最大得分数
                conf.int = F, # 添加置信区间，很难看，可以不要
                funlabel="Risk of DFU")  
plot(nom) 


#还可以进行一些美化，这个函数的参数非常多，大家可以通过帮助文档自己学习：

# plot(nom1,
#      col.grid=c("tomato","grey")
#      #conf.space = c(0.3,0.5) # 置信区间位置
# ) 

target_gene






# 02绘制列线图163-------------------------------
dd_nomo_163 <- dd_163
dd_nomo_163$group <- as.factor(dd_nomo_163$group)

library(rms)

rms_fit_163 <- rms::lrm(group ~ CALCOCO1 + FUCA1 + CCND2 ,
                    data = dd_nomo_163)
#fit1
dd_nomo_163 <- datadist(dd_nomo_163)
options(datadist="dd_nomo_163")
nom_163 <- nomogram(rms_fit_163, fun=plogis,
                fun.at=c(0.001,0.1,0.25,0.5,0.75,0.9,0.99),
                lp=T, # 是否显示线性预测值
                maxscale = 100, # 最大得分数
                conf.int = F, # 添加置信区间，很难看，可以不要
                funlabel="Risk of IPH")  
plot(nom_163) 






