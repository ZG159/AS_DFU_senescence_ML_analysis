expr_134 <- readRDS("../任务：DFU and MS/expr_134.Rdata")
gene <- read.csv("../任务：DFU and MS/SCG1_DFU_IPH_expr134.csv")
gene <- gene$x

gene
#01 lasso---------
library(Ipaper) 
library(glmnet)
library(tibble)
library(dplyr)
set.seed(2025)

## 1.数据的准备及预处理 -------------------------------------------------------------

exp_lasso <- expr_134[gene,]
metadata1 <- data.frame(sample = colnames(exp_lasso), group = c(rep("Con",10),rep("DFU",11)))

exp_lasso <- exp_lasso %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  inner_join(metadata1,by = "sample") %>% 
  column_to_rownames("sample")



x_all <- subset(exp_lasso, select = -group) #基因表达矩阵
y_all <- subset(exp_lasso, select = group)#分组数据


## 2.进行lasso分析 -------------------------------------------------------------

# Lasso筛选变量动态过程图
mod <- glmnet(as.matrix(x_all),as.matrix(y_all), family="binomial",alpha = 1) 
plot(mod,xvar = "lambda", label = T,lwd=1)

#lasso回归最重要的就是选择合适的λ值，可以通过 cv.glmnet函数实现     
cvmod = cv.glmnet(as.matrix(x_all), as.matrix(y_all),family = "binomial")
plot(cvmod, las =1)


# 最终模型的系数估计(使用lambda.minl;也可以使用lambda.1se)
coef.min <-  coef(cvmod, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
#coef.1se <- coef(cvmod, s = "lambda.1se")
active.min <-  which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
#active.1se  <-  which(coef.1se@i != 0) ## 找出那些回归系数没有被惩罚为0的
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
#lasso_geneids <- coef.1se@Dimnames[[1]][coef.1se@i+1] ## 提取基因名称

lasso_geneids <- lasso_geneids[-1] %>% as.data.frame()
colnames(lasso_geneids) <- 'symbol'
lasso_geneids <- lasso_geneids$symbol
lasso_geneids
# lasso结果简单可视化(另一种画法)------------------------
par(mfrow = c(1, 2))
plot(cvmod$glmnet.fit, 
     xvar = 'lambda', 
     label = TRUE, 
     las = 1)
abline(v=log(c(cvmod$lambda.min, cvmod$lambda.1se)), lty="dashed")
plot(cvmod, 
     las =1)

lasso_geneids




