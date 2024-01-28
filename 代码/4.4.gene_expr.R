rm(list=ls())
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)
source("src/functions/v_characteristics_plot_by_group.R")
# 设置结果路径
od <- "results/5.model/checkpoint"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")

checkpoints_data <- readxl::read_xlsx(paste0(Data_Center,"/Gene_Sets/immune_features/checkpiont.xlsx")) %>%
    mutate(Type = case_when(
        Role.with.Immunity %in% c("Activate", "Active") ~ "Activate",
        Role.with.Immunity %in% c("Inhibit") ~ "Inhibit",
        Role.with.Immunity %in% c("TwoSide") ~ "TwoSide"
    )) %>%
    arrange(Type, Symbol)

# 全部的结果
v_characteristics_plot_by_group(
    characteristics_score = train_data$tumor_exprs[rownames(train_data$tumor_exprs) %in% checkpoints_data$Symbol, ] %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("sample"),
    Group =train_model_res$Group,
    od = str_c(od, "/all"), type = "all"
)

# 挑选显著的进行box展示
sign_checkpoints <- read.delim(file = str_glue("{od}/all/SupplementaryTable_all_sign_stat.txt")) %>%
    rownames_to_column("checkpoints") %>%
    filter(pvalue < 0.05) %>%
    pull(checkpoints)
v_characteristics_plot_by_group(
    characteristics_score = train_data$tumor_exprs[rownames(train_data$tumor_exprs) %in% sign_checkpoints, ] %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column("sample"),
    Group =train_model_res$Group,
    od = str_c(od, "/significant_checkpoint"), type = "all"
)
