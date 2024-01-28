rm(list=ls())
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)
source("src/functions/v_characteristics_plot_by_group.R")
# 设置结果路径
od <- "results/9.immune/"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")
# 免疫细胞浸润
data <- readxl::read_xlsx(paste0(Data_Center,'/Gene_Sets/immune_features/28celltype_gene_PMID28052254.xlsx'))
ssgsea_geneSets <- split(as.matrix(data)[, 1], data[, 2])
immunescore <- immune_score(
    exp = train_data$tumor_exprs, arrays = FALSE, perm = 200, group_list = NULL,
    tumor = TRUE, scale_mrna = TRUE, od = od, ssgsea_geneSets = ssgsea_geneSets,
    method = c("ssgsea",'cibersort','estimate')
)
save(immunescore, file = str_glue("{od}/immunescore.RData"))
load(str_glue("{od}/immunescore.RData"))

# Adaptive
type28 <- readxl::read_xlsx(
    str_glue(("{Data_Center}/Gene_Sets/immune_features/28celltype_gene_PMID28052254.xlsx"))
    ,skip=2)
data <- type28 %>% 
    dplyr::filter(Immunity == "Adaptive") %>%
    pull(`Cell type`) %>%
    unique()
dir.create(str_glue("{od}/immune_infiltration"))
conflicts_prefer(clusterProfiler::rename)
conflicts_prefer(base::intersect)
conflicts_prefer(clusterProfiler::select)
v_characteristics_plot_by_group(
    characteristics_score = immunescore$ssgsea %>% select(sample, data),,
    Group = train_model_res$Group,
    od = str_glue("{od}/immune_infiltration/ssgsea_adaptive"),
    type = "ssgsea",
    heatplot_by_scale = TRUE,
    cluster_rows = TRUE
)
# Innate
data <- type28 %>% 
    filter(Immunity == "Innate") %>%
    pull(`Cell type`) %>%
    unique()

v_characteristics_plot_by_group(
    characteristics_score = immunescore$ssgsea %>% select(sample, data),
    Group = train_model_res$Group,
    od = str_glue("{od}immune_infiltration/ssgsea_innate/"),
    type = "ssgsea",
    heatplot_by_scale = TRUE, 
    cluster_rows = TRUE
)
# other
geneSets <- GSEABase::getGmt("/Volumes/T7/DataHub/Gene_Sets/MsigDB/h.all.v2022.1.Hs.symbols.gmt")
hallmarkscore <- immune_score(
    exp = train_data$tumor_exprs, arrays = FALSE, perm = 200, group_list = NULL,
    tumor = TRUE, scale_mrna = TRUE, od = od, ssgsea_geneSets = geneSets,
    method = c("ssgsea")
)
save(hallmarkscore, file = str_glue("{od}/hallmarkscore.RData"))
load(str_glue("{od}/hallmarkscore.RData"))

v_characteristics_plot_by_group(
    characteristics_score = hallmarkscore$ssgsea,
    Group =  train_model_res$Group,
    od = str_glue("{od}hallmark"),
    type = "ssgsea",
    heatplot_by_scale = TRUE,
    cluster_rows = TRUE
)
