rm(list = ls())
source("src/0.config.R")
# 设置结果路径
od <- "results/10.treatment/"
suppressWarnings(dir.create(od, recursive = TRUE))

# 1.读入药物数据
library(readxl)
dat1 <- read_excel(path = "data/DTP NCI-60/DTP_NCI60_ZSCORE.xlsx", skip = 7)

colnames(dat1) <- dat1[1, ]
dat1 <- dat1[-1, -c(67, 68)]

# 筛选药物标准
# 选取经过临床试验（Clinical trial）和FDA批准（FDA approved）的药物结果
dat1 <- dat1[dat1$`FDA status` %in% c("FDA approved", "Clinical trial"), ]
dat1 <- dat1[, -c(1, 3:6)]

fwrite(dat1, paste0(od, "/CellMiner_drug.csv"))

# 3.读入表达数据#
load("data/train_data.RData")
exp <- train_data$tumor_exprs

library(impute)
library(limma)

# 4.整理数据
# 读取药物输入文件
drugDat <- read.table(paste0(od, "/CellMiner_drug.csv"), sep = ",", header = T, check.names = F)
drugDat <- as.matrix(drugDat)
rownames(drugDat) <- drugDat[, 1]
drug <- drugDat[, 2:ncol(drugDat)]
dimnames <- list(rownames(drug), colnames(drug))
data <- matrix(as.numeric(as.matrix(drug)), nrow = nrow(drug), dimnames = dimnames)
head(data)

# 考虑到药物敏感性数据中存在部分NA缺失值，通过impute.knn()函数来评估并补齐药物数据。其中，impute.knn()函数是一个使用最近邻平均来估算缺少的表达式数据的函数。
mat <- impute.knn(data)
drug <- mat$data
drug <- avereps(drug) %>%
    t() %>%
    as.data.frame()

colnames(drug)[1:12]

# 提取特定基因表达
library(WGCNA)
library(tidyr)
inputgene <- c("TP53", "PTEN", "BCAT2", "EGFR", "TMEM178A")
gl <- intersect(inputgene, row.names(exp))
exp <- exp[gl, ] %>%
    t() %>%
    as.data.frame()

identical(rownames(exp), rownames(drug))

dim(drug)
rownames(exp)
### 5.药敏相关性分析
## ======药物敏感性计算
outTab <- data.frame()
for (Gene in row.names(exp)) {
    x <- as.numeric(exp[Gene, ])
    # 对药物循环
    for (Drug in row.names(drug)) {
        y <- as.numeric(drug[Drug, ])
        corT <- cor.test(x, y, method = "spearman")
        cor <- corTSestimate
        pvalue <- corTsp.value
        if (pvalue < 0.01) {
            outVector <- cbind(Gene, Drug, cor, pvalue)
            outTab <- rbind(outTab, outVector)
        }
    }
}

dim(exp)
view(drug) # 60 409
