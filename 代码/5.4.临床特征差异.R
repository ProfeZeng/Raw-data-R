rm(list=ls())
source("")
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)

# 设置结果路径
od <- "results/5.model/clinical"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")
com_gene <- fread("results/3.WGCNA/com_gene2.csv") %>% pull(1)

characteristics_score <- train_data$tumor_exprs[com_gene, ] %>%
    na.omit() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample")

# 不同临床特征分组

source("src/functions/v_characteristics_plot_by_group.R",local = TRUE)
v_characteristics_plot_by_group(
    characteristics_score = characteristics_score,
    feature2show=parameter_list$n_exprs,
    Group = train_model_res$Group,
    od = od # str_glue("{od}/{x}/"),
    , type = "Risk"
)

conflicts_prefer(dplyr::select)
fwrite(as.data.table(group_infor,keep.rownames = "Sample"),str_glue("{od}/cluster.csv"))
fwrite(characteristics_score, str_glue("{od}/exprs.csv.gz"))

