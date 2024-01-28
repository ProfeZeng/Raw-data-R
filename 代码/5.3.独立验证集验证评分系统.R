rm(list=ls())
# 导入函数

source("src/0.config.R")
# 设置结果路径
od = "results/5.model/validate"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("results/5.model/signature/cox_res.RData")
load("data/train_data.RData")

# 验证集
source("src/functions/v_datacenter_validation.R")
signature_coef <- if (str_detect(parameter_list$model, "PCA")) {
    tibble(signature = cox_res$cox_signature, coef = 1)
} else {
    load("results/5.model/signature/signature_coef.RData")
    signature_coef
}
source("src/functions/v_datacenter_validation.R")
DATA_DIR = '/Volumes/T7/DataHub'


valid_model_res <- v_datacenter_validation(
    signature_coef = signature_coef,
    od = od,
    edition="v1",
    cox_tibble = cox_res$coxResult,
    cancer = parameter_list$cancer,
    outcome_list = parameter_list$valid_survival_outcome_list,
    method = parameter_list$model,
    time = c(1,3,5),
    best_cut=F,
    prognostic_independence = T,
    cohort_list=c("LIRI-JP_seq"),
    factors = c("Age" ,"Score","Grade","Stage")
)

save(valid_model_res, file = str_glue("{od}/valid_model_res.RData"))
load("/Volumes/T7/DataHub/ICGC/v1/LIRI-JP_seq/clinical.RData")
colnames(clinical)

grep("UNC5H2",rownames(expression),value = T,ignore.case=T)

rownames(expression) <- str_replace(rownames(expression),"CSPG2","VCAN")
rownames(expression) <- str_replace(rownames(expression),"PSCD3","CYTH3")
rownames(expression) <- str_replace(rownames(expression),"HIST1H4K","H4C12")
save(expression,file="/Volumes/T7/DataHub/GEO/v1/GSE14520_GPL3921/expression.RData")
