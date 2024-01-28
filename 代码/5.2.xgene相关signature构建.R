rm(list = ls())
# 导入函数
source("src/0.config.R")
# 设置结果路径
od <- "results/5.model/signature/"
suppressWarnings(dir.create(od, recursive = TRUE))
conflict_prefer("rename", "dplyr")
# 导入数据
load("data/train_data.RData")
# load("data/xgene.RData")
load("results/5.model/signature/cox_res.RData")
cox_tibble <- cox_res$coxResult
cox_res$cox_signature
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
set.seed(0)
signature_coef <- lasso_model(
    timecol = paste0(parameter_list$train_survival_outcome[1], ".Time"),
    statuscol = paste0(parameter_list$train_survival_outcome[1], ".Status"),
    signaturelist = cox_res$cox_signature,
    od = od,
    seed = 0,
    exp = train_data$tumor_exprs,
    clin = mutate(train_data$data_clinical,sample=Sample),
    signaturetype = "cox_gene"
)
if (str_detect(parameter_list$model, pattern = "PCA")) {
    signature_coef <- tibble(signature = "PC1PC2", coef = 1)
} else {
    exprs <- train_data$tumor_exprs
}
conflict_prefer("ggarrange", "ggpubr")
source('src/functions/v_modelscore_km_roc.R')
train_model_res <- modelscore_km_roc(
    timecol = paste0(parameter_list$train_survival_outcome[1], ".Time"),
    statuscol = paste0(parameter_list$train_survival_outcome[1], ".Status"),
    signature_coef = signature_coef,
    od = str_glue("{od}/train/"),
    # roc_time = parameter_list$roc_time,
    no_roc_pheatmap = F,
    exp = exprs,
    best_cut=F,
    clin = mutate(train_data$data_clinical,sample=Sample),
    dataset = parameter_list$train_cohort,
    xlab = str_c(parameter_list$train_survival_outcome, " days")
)

save(train_model_res, file = str_glue("{od}/train_model_res.RData"))
save(signature_coef, file = str_glue("{od}/signature_coef.RData"))
# 评分系统的预后独立性
source("src/functions/v_uni_multi_cox.R")
res <- uni_multi_cox(
    timecol = paste0(parameter_list$train_survival_outcome[1], ".Time"),
    statuscol = paste0(parameter_list$train_survival_outcome[1], ".Status"),
    od = str_glue("{od}/train/"),
    infor = train_model_res$Group %>% rownames_to_column("sample") %>% merge(., train_data$data_clinical) %>% rename(Score = Group),
    factors = c("Sex" ,"Age" ,"Score"),
    dataset = parameter_list$train_cohort, coxp = 1, w = 10
)
colnames(train_data$data_clinical)

# score = as_tibble(train_model_res$Group,rownames = "sample") %>% 
#     rename(Score = Group) 
# nrow(score)
# expr = train_data$tumor_exprs[signature_coef$signature,] %>% 
#     t() %>% 
#     as_tibble(.,rownames = "sample") %>% 
#     mutate(INCENP=ifelse(INCENP>median(INCENP),"High","Low")
#         ,TP53=ifelse(TP53>median(TP53),"High","Low")
#         ,NUP153=ifelse(NUP153>median(NUP153),"High","Low")
#         ,UBA2=ifelse(UBA2>median(UBA2),"High","Low")
#         ,NUP43=ifelse(NUP43>median(NUP43),"High","Low")
#     ) %>% 
#     mutate(INCENP=factor(INCENP,levels = c("Low","High"))
#     ,TP53=factor(TP53,levels = c("Low","High"))
#     ,NUP153=factor(NUP153,levels = c("Low","High"))
#     ,UBA2=factor(UBA2,levels = c("Low","High"))
#     ,NUP43=factor(NUP43,levels = c("Low","High"))
#     )
# cc=train_data$data_clinical %>% 
#     dplyr::select(Sample,OS.Time,OS.Status) %>% 
#     rename(sample = Sample) 

# infor=reduce(list(score,expr,cc),inner_join,by='sample')
# head(infor)

# source("src/functions/v_uni_multi_cox.R")
# res <- v_uni_multi_cox(
#     timecol = paste0(parameter_list$train_survival_outcome[1], ".Time"),
#     statuscol = paste0(parameter_list$train_survival_outcome[1], ".Status"),
#     od = str_glue("{od}/train/"),
#     infor = infor,
#     factors = c("Score",signature_coef$signature),
#     dataset = parameter_list$train_cohort, coxp = 0.05, w = 10
# )
