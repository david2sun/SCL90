library(stringr)
df <- read.csv("/Users/sunyajun/Desktop/SCL-90分析/scl-90_1.csv",
               header = T, sep = ",", fileEncoding = "gb2312", stringsAsFactors = F)
names(df)[18:107] <- paste0("y",1:90)
##参考标准：总分超过 160 分  ###################
##########  或阳性项目数超过 43 项，############
##########  或任一因子分超过 2 分，可考虑筛查阳性，须进一步检查。###
# 总均分 raw Global Severity Index (GSI)
df$total_scores <- apply(as.matrix(df[,18:107]), MARGIN = 1, sum) / 90
c(mean(df$total_scores), sd(df$total_scores)) # 均数和标准差
# convert raw GSI into T-scores
# T = (Z x 10) + 50
df$GSI_Tscore <- scale(df$total_scores) * 10 + 50
mean(scale(df$total_scores))
# sum(df$GSI_Tscore >= 63)
# mean(df$GSI_Tscore >= 63)
# library(samplingbook)
# Sprop(df$GSI_Tscore >= 63, N=23690) ## 样本率估计
# # df$group[df$GSI_Tscore >= 63] <- "high risk"
# # df$group[df$GSI_Tscore < 63] <- "low risk"
# df$group[df$GSI_Tscore >= 63] <- 1
# df$group[df$GSI_Tscore < 63] <- 0
# 阳性项目数 Positive Symptom Total (PST)
df$pos_items <- apply(as.matrix(df[,18:107]) > 1, MARGIN = 1, sum)
c(mean = mean(df$pos_items[df$pos_items > 0]), sd = sd(df$pos_items[df$pos_items > 0]), median = median(df$pos_items[df$pos_items > 0])) 
# 阳性症状均分 Positive Symptom Distress Index (PSDI)
z <- as.matrix(df[,18:107])
df$pos_score <- apply(ifelse(z > 1, z, 0), MARGIN = 1, sum) # 计算阳性分数
pos_mean <- df$pos_score[df$pos_score > 1] / df$pos_items[df$pos_score > 1]
c(mean = mean(pos_mean), sd = sd(pos_mean), median = median(pos_mean)) 
# 躯体因子
quti <- c(1, 4, 12, 27, 40, 42, 48, 49, 52, 53, 56, 58)
df$factor_quti <- apply(as.matrix(df[,18:107])[,quti], MARGIN = 1, sum)
df$factor_quti_score <- df$factor_quti / length(quti)
# 强迫因子
qiangpo <- c(3,9,10,28,38,45,46,51,55,65)
df$factor_qiangpo <- apply(as.matrix(df[,18:107])[,qiangpo], MARGIN = 1, sum)
df$factor_qiangpo_score <- df$factor_qiangpo / length(qiangpo)
# 人际关系敏感
rela_se <- c(6,21,34,36,37,41,61,69,73)
df$factor_rela_se <- apply(as.matrix(df[,18:107])[,rela_se], MARGIN = 1, sum)
df$factor_rela_se_score <- df$factor_rela_se / length(rela_se)
# 抑郁
depression <- c(5,14,15,20,22,26,29,30,31,32,54,71,79)
df$factor_depression <- apply(as.matrix(df[,18:107])[,depression], MARGIN = 1, sum)
df$factor_depression_score <- df$factor_depression / length(depression)
# 焦虑
anxiety <- c(2,17,23,33,39,57,72,78,80,86)
df$factor_anxiety <- apply(as.matrix(df[,18:107])[,anxiety], MARGIN = 1, sum)
df$factor_anxiety_score <- df$factor_anxiety / length(anxiety)
# 敌对
hostility <-  c(11,24,63,67,74,81)
df$factor_hostility <- apply(as.matrix(df[,18:107])[,hostility], MARGIN = 1, sum)
df$factor_hostility_score <- df$factor_hostility / length(hostility)
# 恐怖
terror <- c(13,25,47,50,70,75,82)
df$factor_terror <- apply(as.matrix(df[,18:107])[,terror], MARGIN = 1, sum)
df$factor_terror_score <- df$factor_terror / length(terror)
# 偏执
paranoid <- c(8,18,43,68,76,83)
df$factor_paranoid <- apply(as.matrix(df[,18:107])[,paranoid], MARGIN = 1, sum)
df$factor_paranoid_score <- df$factor_paranoid / length(paranoid)
# 精神病性
Psyc <- c(7,16,35,62,77,84,85,87,88,90)
df$factor_Psyc <- apply(as.matrix(df[,18:107])[,Psyc], MARGIN = 1, sum)
df$factor_Psyc_score <- df$factor_Psyc / length(Psyc)
# 其他
other <- c(quti,qiangpo,rela_se,depression,anxiety,hostility,
           terror,paranoid,Psyc)
df$factor_other <- apply(as.matrix(df[,18:107])[,-other], MARGIN = 1, sum)
df$factor_other_score <- df$factor_other / (90 - length(other))
########################## 筛查阳性判断 ########################
factors <- grep("factor.*score", names(df), fixed = FALSE, value = T)
df$res <- with(data = df,
          ifelse(total_scores * 90 > 160  | pos_items > 43 | apply(df[factors] > 2, 1, any), 1, 0))
sum(df$res == 1) ## 阳性总数
mean(df$res == 1) ## 阳性比例
# write.csv(df, file = "final.csv", row.names = F, fileEncoding = "GB2312")

# 计算10个因子分 均值（标准值）x ± s
# factors <- grep("factor.*score", names(df), fixed = FALSE, value = T) # regular expression
mydf <- df[factors]
sapply(mydf, function(x) c(mean = mean(x),sd = sd(x)))
sapply(mydf, function(x) c(n = sum(x <= 2), prop = mean(x <= 2)))
# 总均分与常模t检验（单样本）
t.test(df$total_scores, mu = 1.44)
t.test(log(df$total_scores), mu = log(1.44))
qqPlot(log(df$total_scores))
# 9个因子分与常模t检验（单样本）
mu <- c(1.37,1.62,1.65,1.50,1.39,1.46,1.23,1.43,1.29)
mylist <- list()
for (i in 1:9) {
  mylist[[i]] <- t.test(df[,factors[i]], mu = mu[i])
}
mylist
for (i in 1:9) {
  mylist[[i]] <- t.test(log(df[,factors[i]]), mu = log(mu[i]))
}
mylist
# nine dimension severity distrubution 
# library(gmodels)
res1 <- list()
for (i in 1:length(factors)) {
  var <- cut(df[, factors[i]], breaks = c(0,1,2,3,4,5))
  res1[[factors[i]]] <- rbind(n = table(var), prop = round(100*prop.table(table(var)), 2))
  }
res1
# 年龄分组
# df$age.group <- cut(df$X4.您的年龄., breaks = c(0, 20, 30, 40, 50, 60, 100),
#                     right = F)
df$age.group <- cut(df$X4.您的年龄., breaks = c(0,30,40,50,100),
                    right = F)
levels(df$age.group)
table(df$age.group)
round(100 * prop.table(table(df$age.group)),2) # 计算比例
# 不同年龄组筛查阳性率比较
addmargins(table(df$age.group, df$res))
prop.table(table(df$age.group, df$res), margin = 1)
chisq.test(table(df$age.group, df$res))
library(car);library(gmodels);library(rcompanion)
chisq.test(table(df$age.group, df$res))
summary(table(df$age.group, df$res))
y1 <- pairwiseNominalIndependence(table(df$age.group, df$res), method = "bonferroni")
cldList(comparison = y1$Comparison, 
        p.value    = y1$p.adj.Fisher, 
        threshold  = 0.05)
## 不同性别阳性率比较 ###
tbl.sex <- table(df$X3.您的性别., df$res)
addmargins(tbl.sex)
round(100 * prop.table(tbl.sex, margin = 1), 2)
chisq.test(tbl.sex)
## 不同婚姻状态阳性率比较 ###
df$X5.您的婚姻状况. <- ifelse(df$X5.您的婚姻状况. == 2, df$X5.您的婚姻状况., 1) 
tbl.mari <- table(df$X5.您的婚姻状况., df$res)
addmargins(tbl.mari)
round(100 * prop.table(tbl.mari, margin = 1), 2)
chisq.test(tbl.mari)
## 不同受教育程度阳性率比较
df$X6.您的受教育程度.[df$X6.您的受教育程度. == 1] <- 2 # 合并高中
tbl <- table(df$X6.您的受教育程度., df$res)
addmargins(tbl)
round(100 * prop.table(tbl, margin = 1), 2)
chisq.test(tbl)
y1 <- pairwiseNominalIndependence(tbl, method = "bonferroni")
cldList(comparison = y1$Comparison, 
        p.value    = y1$p.adj.Fisher, 
        threshold  = 0.05)
## 不同工作年限阳性率比较
df$work.years <- cut(df$X8.您的工作年限.年..不满1年按1年算., 
                     breaks = c(0,5,10,20,60), right = F)
tbl <- table(df$work.years, df$res)
addmargins(tbl)
round(100 * prop.table(tbl, margin = 1), 2)
chisq.test(tbl)
y1 <- pairwiseNominalIndependence(tbl, method = "bonferroni")
cldList(comparison = y1$Comparison, 
        p.value    = y1$p.adj.Fisher, 
        threshold  = 0.05)
## 不同工作经历阳性率比较 ###
tbl <- table(df$X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等., df$res)
addmargins(tbl)
round(100 * prop.table(tbl, margin = 1), 2)
chisq.test(tbl)

# 不同年龄组GSI T-scores比较
# table(df$age.group)
# library(reshape)
# stasfun <- function(x) c(n = length(x), mean = mean(x), sd = sd(x), median = median(x))
# df %>% melt(., measure.vars = c("GSI_Tscore")) %>%
#        cast(., age.group ~ ., stasfun) %>%
#        dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2)) # 分组计算均数和标准差
# df %>% melt(., measure.vars = c("total_scores")) %>%
#   cast(., age.group ~ ., stasfun) %>%
#   dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2))
# # 方差分析
# library(car)
# qqPlot(GSI_Tscore ~ age.group, data = df) #正态性评估
# # qqPlot(lm(GSI_Tscore ~ age.group, data = df)) 方差分析的正态性假设检验
# library(nortest)
# lillie.test(df$GSI_Tscore)  # (Kolmogorov-Smirnov) normality test
# # 基于分组的正态性检验
# nor.test <- function(x) {
#   result <- shapiro.test(x)
#   c(W = round(result$statistic, 4), p = round(result$p.value, 4))
# }
# aggregate(GSI_Tscore ~ age.group, data = df, FUN = nor.test)
# 
# 
# 
# # 提示不具有正态性
# car::leveneTest(lm(GSI_Tscore ~ age.group, data = df)) #方差齐性检验
# # 提示方差不齐
# # Welch ANOVA（方差不齐时使用）
# # oneway.test(GSI_Tscore ~ age.group, data = df)
# # 多组的非参数检验（成组设计的秩和检验）
# kruskal.test(GSI_Tscore ~ age.group, data = df)
# # 多组两两比较
# source("http://www.statmethods.net/RiA/wmc.txt")
# wmc(GSI_Tscore ~ age.group, data = df, method = "holm")
# # 不同性别GSI T-scores比较
# df %>% melt(., measure.vars = c("GSI_Tscore")) %>%
#   cast(., X3.您的性别. ~ ., stasfun) %>%
#   dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2)) # 分组计算均数和标准差
# prop.table(table(df$X3.您的性别.)) # 比例
# qqPlot(GSI_Tscore ~ X3.您的性别., data = df) # 正态性评估
# ## 结果提示非正态性
# wilcox.test(GSI_Tscore ~ X3.您的性别., data = df) # stats包函数
# coin::wilcox_test(GSI_Tscore ~ factor(X3.您的性别.), data = df)
# # W = 2631874, p-value < 2.2e-16
# # 不同婚姻状况GSI T-scores比较
# df$X5.您的婚姻状况. <- ifelse(df$X5.您的婚姻状况. == 2, df$X5.您的婚姻状况., 1) # 两分类
# df %>% melt(., measure.vars = c("GSI_Tscore")) %>%
#   cast(., X5.您的婚姻状况. ~ ., stasfun) %>%
#   dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2)) # 分组计算均数和标准差
# round(100*prop.table(table(df$X5.您的婚姻状况.)), 2)
# coin::wilcox_test(GSI_Tscore ~ factor(X5.您的婚姻状况.), data = df)
# wilcox.test(GSI_Tscore ~ X5.您的婚姻状况., data = df) # 两样本
# kruskal.test(GSI_Tscore ~ X5.您的婚姻状况., data = df) # 多样本总体检验
# # Kruskal-Wallis chi-squared = 3.8845, df = 3, p-value = 0.2742
# wmc(GSI_Tscore ~ X5.您的婚姻状况., data = df, method = "holm") ##两两比较
# # 不同受教育程度GSI T-scores比较
# df$X6.您的受教育程度.[df$X6.您的受教育程度. == 1] <- 2 # 合并高中
# df %>% melt(., measure.vars = c("GSI_Tscore")) %>%
#   cast(., X6.您的受教育程度. ~ ., stasfun) %>%
#   dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2)) # 分组计算均数和标准差
# round(100*prop.table(table(df$X6.您的受教育程度.)), 2)
# kruskal.test(GSI_Tscore ~ X6.您的受教育程度., data = df) 
# # Kruskal-Wallis chi-squared = 13.73, df = 3, p-value = 0.003297
# wmc(GSI_Tscore ~ X6.您的受教育程度., data = df, method = "holm")
# # 不同工作年限GSI T-score比较
# df$work.years <- cut(df$X8.您的工作年限.年..不满1年按1年算., 
#                      breaks = c(0,5,10,20,60), right = F)
# df %>% melt(., measure.vars = c("GSI_Tscore")) %>%
#   cast(., work.years ~ ., stasfun) %>%
#   dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2)) # 分组计算均数和标准差
# round(100*prop.table(table(df$work.years)), 2)
# kruskal.test(GSI_Tscore ~ work.years, data = df) 
# # Kruskal-Wallis chi-squared = 16.373, df = 3, p-value =
# #   0.0009508
# wmc(GSI_Tscore ~ work.years, data = df, method = "holm") ##两两比较
# # 不同COVID-19工作经历GSI T-score比较
# df %>% melt(., measure.vars = c("GSI_Tscore")) %>%
#   cast(., X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等. ~ ., stasfun) %>%
#   dplyr::mutate(mean = round(mean, 2), sd = round(sd, 2), median = round(median, 2)) # 分组计算均数和标准差
# round(100*prop.table(table(df$X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等.)), 2)
# wilcox.test(GSI_Tscore ~ X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等., data = df)
# coin::wilcox_test(GSI_Tscore ~ factor(X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等.), data = df)
#高危发生影响因素分析（单因素）
table(df$age.group, df$group) # 年龄组卡方检验
addmargins(table(df$age.group, df$group))
round(100*prop.table(table(df$age.group, df$group), margin = 1), 2)
coin::chisq_test(table(df$age.group, df$group))
table(df$X3.您的性别., df$group) # 性别组卡方检验
addmargins(table(df$X3.您的性别., df$group))
round(100*prop.table(table(df$X3.您的性别., df$group), margin = 1), 2)
coin::chisq_test(table(df$X3.您的性别., df$group))
table(df$X5.您的婚姻状况., df$group) # 婚姻组卡方检验
addmargins(table(df$X5.您的婚姻状况., df$group))
round(100*prop.table(table(df$X5.您的婚姻状况., df$group), margin = 1), 2)
coin::chisq_test(table(df$X5.您的婚姻状况., df$group))
table(df$X6.您的受教育程度., df$group) # 教育程度卡方检验
addmargins(table(df$X6.您的受教育程度., df$group))
round(100*prop.table(table(df$X6.您的受教育程度., df$group), margin = 1), 2)
coin::chisq_test(table(df$X6.您的受教育程度., df$group))
table(df$work.years, df$group) # 工作年限卡方检验
addmargins(table(df$work.years, df$group))
round(100*prop.table(table(df$work.years, df$group), margin = 1), 2)
coin::chisq_test(table(df$work.years, df$group))
table(df$X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等., df$group) # COVID-19工作经历卡方检验
addmargins(table(df$X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等., df$group))
round(100*prop.table(table(df$X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等., df$group), margin = 1), 2)
coin::chisq_test(table(df$X9.您是否从事新冠肺炎疫情相关工作..如与确诊.疑似患者.密切接触者.隔离人员.流动人员排查相关工作等., df$group))
# logistic回归多因素分析
# names(df)[c(10:13,15,16)] <- c("sex","age","marriage","education","occup.years","covid.exp")
# # df$covid.exp <- ifelse(df$covid.exp == 2, 0, df$covid.exp)
# df$covid.exp <- relevel(factor(df$covid.exp), ref = "2")
# df$sex <- relevel(factor(df$sex), ref = "1")
# # df$sex <- ifelse(df$sex == 1, df$sex, 0)
# # df$sex <- factor(df$sex, levels = c("1","0"), labels = c("男","女"))
# df$marriage <- factor(df$marriage)
# df$education <- factor(df$education)
# df$group <- factor(df$group)
# levels(df$age.group)
# df$age.group <- relevel(df$age.group, ref = "[50,100)")
# df$work.years <- cut(df$occup.years, breaks = c(0,5,10,20,60), right = F)
# levels(df$work.years)
# df$work.years <- relevel(df$work.years, ref = "[0,5)")
# ## 看看age和occup.years相关性
# cor.test(df$age, df$occup.years, method = "spearman") # 结果提示有高度相关性
# fit1 <- glm(group ~ sex + age.group + marriage + education + covid.exp, 
#            data = df, family = binomial())
# summary(fit1)
# exp(coef(fit1))
# fit2 <- glm(group ~ sex + age.group + marriage + education + work.years + covid.exp, 
#             data = df, family = binomial())
# fit3 <- glm(group ~ sex + marriage + education  + covid.exp, 
#             data = df, family = binomial())
## logistic回归2 以阳性率为结局
names(df)[c(10:13,15,16)] <- c("sex","age","marriage","education","occup.years","covid.exp")
# df$covid.exp <- ifelse(df$covid.exp == 2, 0, df$covid.exp)
df$covid.exp <- relevel(factor(df$covid.exp), ref = "2")
df$sex <- relevel(factor(df$sex), ref = "1")
# df$sex <- ifelse(df$sex == 1, df$sex, 0)
# df$sex <- factor(df$sex, levels = c("1","0"), labels = c("男","女"))
df$marriage <- factor(df$marriage)
df$education <- factor(df$education)
df$res <- factor(df$res)
levels(df$age.group)
df$age.group <- relevel(df$age.group, ref = "[50,100)")
df$work.years <- cut(df$occup.years, breaks = c(0,5,10,20,60), right = F)
levels(df$work.years)
df$work.years <- relevel(df$work.years, ref = "[0,5)")
cor.test(df$age, df$occup.years, method = "spearman") # 结果提示有高度相关性
fit1 <- glm(res ~ sex + age.group + marriage + education + work.years + covid.exp, 
            data = df, family = binomial())
round(cbind(summary(fit1)$coefficients, 
      coef = exp(coef(fit1)), 
      exp(confint(fit1))), 4)
fit2 <- step(fit1) ## 去掉无意义的变量
summary(fit2)
round(cbind(coef = exp(coef(fit2)), exp(confint(fit2))), 4) ## 点估计及95%CI
anova(fit2, test = "Chisq") ##relative importance analysis based on residual.deviance
library(dominanceanalysis) 
dapres<-dominanceAnalysis(fit2)
averageContribution(dapres,fit.functions = "r2.m")
## 年龄组和工作年限 重新分类 ###
df$age.group2 <- ifelse(df$age < 50, 1, 0)
df$work.years2 <- ifelse(df$occup.years >= 5, 1, 0)
fit3 <- glm(res ~ sex + age.group2 + marriage + work.years2 + covid.exp, 
            data = df, family = binomial())
round(cbind(summary(fit3)$coefficients, 
            coef = exp(coef(fit3)), 
            exp(confint(fit3))), 4)
dapres<-dominanceAnalysis(fit3)
x <- averageContribution(dapres,fit.functions = "r2.m")
round(sort(x$r2.m), 4)
anova(fit3, test = "Chisq")
dominanceMatrix(dapres, type="complete",fit.functions = "r2.m", ordered=TRUE)
dominanceMatrix(dapres, type="conditional",fit.functions = "r2.m", ordered=TRUE)
dominanceMatrix(dapres, type="general",fit.functions = "r2.m", ordered=TRUE)
plot(dapres, which.graph ="general",fit.function = "r2.m") 
library(ggplot2)
df.gg <- data.frame(name=c("女性", "年龄<50岁", "已婚", "工作年限≥5年" ,"新冠疫情工作经历"), value = x$r2.m)
ggplot(data = df.gg, aes(x=reorder(name, value, dplyr::desc), y=value)) +
      geom_bar(aes(x=reorder(name, value, dplyr::desc), y=value, fill=name), stat = "identity") + 
      labs(x="变量", y=expression(bold(R[McFadden]^{'2'})), 
           fill="变量"#, caption.title = "医务人员心理异常筛查阳性Logistic回归模型变量相对重要性") + 
      ) +
      geom_text(stat="identity",aes(label=round(value, 4)), 
                  vjust = 0, family = "STKaiti") +
      scale_fill_discrete(breaks=levels(reorder(df.gg$name, df.gg$value, dplyr::desc))) +
  theme_bw(base_size = 12, base_family = "STKaiti") +
  theme(#plot.caption  = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text  =  element_text(face="bold", size = 11),
            axis.title = element_text(face = "bold", size = 12)) 
       
  
      


