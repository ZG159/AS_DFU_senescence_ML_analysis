set.seed(2025)   # 为了保持结果的一致性
library(e1071)
library(caret)#集成了上百种分类和回归算法
library(kernlab)
library(ggplot2)
#1.前期的数据整理-------------
exp_svm <- exp_lasso


x_all <- subset(exp_svm, select = -group) #基因表达矩阵
y_all <- subset(exp_svm, select = group)#分组数据
y_all$group <- as.factor(y_all$group)


# 2.设置重复交叉验证------------
control <- trainControl(method = "repeatedcv",   
                        # 交叉验证方法,这里选择重复的k折交叉验证
                        number = 10,     
                        # k值,将数据分成5份,可以根据样本量选择5-10
                        repeats = 5,     
                        # 重复次数,可选2-5次,这里选择1次
                        search = "random"   
                        # 选择重复交叉验证时,每轮选择 folds 的方式,可选"systematic"连续选择与"random"随机选择,这里选择随机     
)


#3.使用 SVM-RFE 特征选择方法---------------


svm_rfe <- rfe(x_all,#输入待筛选变量的表达矩阵
               y_all$group,#输入因变量
               #,sizes参数用于指定我们希望选择的特征数量。
               #它决定了RFE算法会考虑选择哪些特征数量的方案。
               #比如sizes = c(1:5, 10)，则RFE算法会考虑选择1个特征,2个特征,3个特征,4个特征,5个特征和10个特征的特征选择方案。
               sizes = 1:15,
               #选择1-10个特征
               rfeControl = rfeControl(functions = caretFuncs,#使用caret包中的评价指标
                                       method = "repeatedcv",# 交叉验证方法,这里选择重复的k折交叉验证
                                       number = 5,          #可以根据样本量选择5-10
                                       repeats = 1, #重复次数，推荐为5，重复5次
                                       verbose = FALSE),
               method = "svmLinear", #使用SVM分类器
               trControl = control,
               preProc = c("center", "scale") #特征预处理,标准化
               
)
importance <- varImp(svm_rfe, scale=FALSE)  #得到基因的权重性，为了后面图形的绘制
importance
svm_rfe$results#可以看到在不同Variables下的Accuracy 
#svm_rfe$results$Variables==8时，准确度最高！
svm_geneids <- rownames(importance)[1:8]


#4.图形绘制--------------- 
plot(svm_rfe,type=c("g", "o"),xlim = c(1,20),ylim = c(0.65,0.80))

predictors(svm_rfe)


svm_geneids
lasso_geneids
rf_geneids

target_gene <- intersect(svm_geneids,intersect(lasso_geneids,rf_geneids))
target_gene



