rm(list=ls())
# 导入函数
source("src/0.config.R")
# 设置结果路径
od = "results/2.DEG"
suppressWarnings(dir.create(od,recursive=TRUE))
# 导入数据
load("data/train_data.RData")
load('data/xgene.RData')

train_data$data_pdata$Group %>% table()

train_data$data_exprs %>% ncol() # 419
train_data$tumor_exprs %>% ncol() # 363
length(xgene) # 85
deg_res <- limma_deg(
    od = od,
    DEG_exp = train_data$data_exprs, 
    DEG_pdata = train_data$data_pdata,
    controlLabel = "Normal",
    caseLabel = 'Tumor',
    DEG_FC = parameter_list$log2fc, 
    DEG_P = parameter_list$deg_p, 
    pvalue = NULL, 
    saveplot = T, 
    color_fun = color_fun1
)
# Tumor_vs_Normal
# log2FC=0.585,FDR=0.05
# DEGs up down
# 8831 7793 1038
degs <- deg_res$DEGs
library(ggvenn)
ggvenn(list(DEGs=deg_res$DEG,CRRGs=xgene),show_percentage=F)
ggsave(paste0(od,'/venn.pdf'),width = 6,height = 6)

com_gene=intersect(deg_res$DEG,xgene)
save(com_gene,file=paste0(od,'/com_gene.RData'))
