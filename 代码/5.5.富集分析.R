rm(list=ls())
source("")
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)

# 设置结果路径
od <- "results/5.model/DEG"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")
com_gene <- fread("results/3.WGCNA/com_gene2.csv") %>% pull(1)

# 
deg_res <- limma_deg(
        od = od, DEG_exp = train_data$tumor_exprs,
        DEG_pdata = rownames_to_column(train_model_res$Group,"sample"),
        controlLabel = "Low", caseLabel = "High",
        DEG_FC = parameter_list$log2fc, DEG_P = parameter_list$deg_p, 
        pvalue = NULL, saveplot = FALSE, color_fun = color_fun1
    )

# High_vs_Low
# log2FC=0.585,FDR=0.05
# DEGs up down
# 2720 2232 488

save(deg_res, file = str_glue("{od}/deg_gene.RData"))

conflicts_prefer(base::setdiff)
up <- deg_res$nrDEG %>% dplyr::filter(Diff=="up") %>% rownames()
fun_res <- enrich(genetype = "up_degs", genelist = up, od = od)

down <- deg_res$nrDEG %>% dplyr::filter(Diff=="down") %>% rownames()
fun_res <- enrich(genetype = "down_degs", genelist = down, od = od)
