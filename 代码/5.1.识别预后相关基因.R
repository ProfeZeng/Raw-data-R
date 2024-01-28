rm(list=ls())
# 导入函数
source("src/0.config.R")
# 设置结果路径
od = "results/5.model/signature"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
xgene <- fread("results/3.WGCNA/com_gene2.csv") %>% pull(1)
# load('results/2.model/deg/deg_gene.RData')
# 单因素cox

cox_res <- signature_cox(
    timecol = paste0(parameter_list$train_survival_outcome[1],".Time"),
    statuscol = paste0(parameter_list$train_survival_outcome[1],".Status"),
    signaturelist = xgene,
    exp = train_data$tumor_exprs, 
    clin = mutate(train_data$data_clinical,sample=Sample), 
    coxp = 0.05, 
    bygroup =F, 
    xlab = str_c(parameter_list$train_survival_outcome, " days"), 
    savekmplot = FALSE
)
# 40
# top20预后基因森林图
# 森林图
forest_plot(
    od = od, input = cox_res$sigcoxResult %>% rownames_to_column("gene") %>% arrange(pvalue), plotN = 20, 
    dataset = parameter_list$train_cohort, h = 8, w = 8,
    signaturecol = "gene", pvaluecol = "pvalue", HRcol = "HR", lower95col = "Low 95%CI",
    upper95col = "High 95%CI"
)
write.table(cox_res$sigcoxResult, file = str_glue("{od}/SupplementaryTable_cox_signature_coxResult.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
save(cox_res, file = str_glue("{od}/cox_res.RData"))

# top6预后基因生存曲线
# length(cox_res$cox_signature)
# row_n <- 2
# col_n <- 3
# p <- survminer::arrange_ggsurvplots(cox_res$kmplot[1:(row_n * col_n)], newPage = F, nrow = row_n, ncol = col_n)
# plotout(p = p, od = od, name = "deg_gene_cox_kmplot", w = col_n * 5, h = row_n * 5)
# dev.off()




