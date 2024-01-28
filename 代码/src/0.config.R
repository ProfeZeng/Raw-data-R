Code_Center = '/Volumes/T7/CodeHub/RCode'
Data_Center = '/Volumes/T7/DataHub'
libSources <- list.files(Code_Center, recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
invisible(lapply(libSources,source))

# 本项目所有提前载入的包和设置
options(bitmapType = "cairo") # for VS code plot
# pre-library packages
library(tidyverse)
library(data.table)
library(magrittr)
library(crayon)
library(conflicted)
conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)


parameter_list <- list(
  #----------①数据整理、SNV、CNV展示、模型效能相关参数----------#
  #-----以下参数需要修改
  cancer = c("LIHC"), # 分析使用的癌症，可以给单个字符也可以给字符向量；当为字符向量时，第一个为训练集所使用的癌症类型，全部为验证集所使用的癌症类型
  version = "v1", # idata版本信息，可选v1或v2，只涉及到验证
  clinical_feature_list = c("Sex",'Stage','Age'), # 用于比较基因表达差异和独立预后验证的临床分组信息
  n_exprs = 30, # 用于画图的基因数目，可以写具体数字数值或"all"，根据xgene的数目自行判断
  n_gene_cell = 30, # 用于画图的基因数目，可以写具体数字数值或"all"，根据xgene的数目自行判断
  n_snv = 30, # 用于画图的基因数目，可以写具体数字数值或"all"，根据xgene的数目自行判断
  n_cnv = 30, # 用于画图的基因数目，可以写具体数字数值或"all"，根据xgene的数目自行判断
  #-----以下参数一般不需要修改
  train_cohort = "TCGA", # 字符，训练集数据库
  train_survival_outcome = "OS", # 字符，训练集使用的临床随访信息终点，只能有一个
  valid_survival_outcome_list = c("OS"), # 字符向量，验证集使用的临床随访信息终点，可以有一个也可以有多个
  drug_list = NULL, # 字符向量，用于化疗耐药性分析的药物名，默认NULL；如果指定的话，则只对指定药物计算和绘图
  cohort_list = NULL, # 字符向量，验证集队列向量，默认NULL；如果指定的话，则只对指定队列进行验证

  #----------②第一步测试precluster参数，一般都不需要修改----------#
  seed = 1314, # 随机数种子
  maxK = 5, # 最大聚类数
  input_distance = c("euclidean", "pearson", "spearman","maximum"), # 聚类的距离计算方式，用于循环,"binary"
  input_clusterAlg = c("hc", "km", "pam"), # 聚类方法，用于循环
  xgene_select_k = c("k2", "k3", "k4"), # 选择的聚类数，筛选聚类结果分类数目

  #----------③第二步测试premodel参数----------#
  #-----以下参数需要修改
  cluster_row_number = c(1), # precluster以后结果文件filter_xgene_cc_res.txt中对应的cluster_row_number列，根据聚类结果图片进行选择
  model = "lasso", # 建模方法，可以选择PCA1、PCA2、lasso、multicox
  gene_cox_bygroup = TRUE, # cox回归时，基因表达量是否要分高低两组
  #-----以下参数一般不需要修改
  logfc_list = c(.585, 1), # 差异分析时，log2FC的取值，用于循环
  degp_list = c(.05, .01), # 差异分析时，p的取值，用于循环
  coxp_list = c(.05, .01), # 单因素cox回归时，p的取值，用于循环
  roc_time = c(1, 3, 5), # 预后模型中的ROC时间1，3，5年

  #----------④根据测试结果选择展示最终结果所需参数----------#
  select_k = "k2", # 聚类数
  select_clusterAlg = "km", # 聚类方法
  select_distance = "euclidean", # 聚类的距离计算方式
  log2fc = 0.585, # 差异分析时log2FC的取值
  deg_p = 0.05, # 差异分析时p的取值
  coxp = .01 # 单因 素cox回归时p的取值
)
color_fun1 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
