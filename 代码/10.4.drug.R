rm(list = ls())
source("src/0.config.R")
# 设置结果路径
od <- "results/10.treatment/"
suppressWarnings(dir.create(od, recursive = TRUE))

# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")
conflict_prefer("between", "dplyr", quiet = TRUE)


install.packages("oncoPredict")
BiocManager::install("GenomicFeatures")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("sva")

data_ic50 <- oncoPredict_calcPhenotype(
    valid_exp = train_data$tumor_exprs, 
    database = "GDSC2", drug_list = NULL, train_exp = NULL, 
    train_ptype = NULL, 
    od = '/Users/victor'
    )

IC50 <- fread(paste0(od,"/DrugPredictions.csv")) %>% 
    as.data.frame() %>%
    mutate_if(is.numeric, log2) %>%
    rename(sample = 1) %>%
    mutate(sample = substr(sample, 1, 16)) %>%
    column_to_rownames("sample")
names(IC50)
library(psych)
samples <- intersect(rownames(IC50), rownames(train_model_res$Group))
IC50 <- as_tibble(IC50[samples, ], rownames = "sample")
IC50[1:3,1:3]
source("src/functions/v_characteristics_plot_by_group.R")
head(IC50[,1:5])

v_characteristics_plot_by_group(
    characteristics_score = IC50,
    Group = train_model_res$Group, od = od, type = "GDSC2",feature2show=6,
)

