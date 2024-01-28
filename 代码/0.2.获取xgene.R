# 设置结果路径
od <- "data"
# -------------------------------获取xgene-------------------------------
# *******目的：从下载的基因集中，得到基因的symbol向量（存为xgene）*******

# 以下是参考代码 分别是https://www.gsea-msigdb.org/gsea/index.jsp
# 和http://amigo.geneontology.org/amigo 数据库下载的数据处理
# 1.msigdb
library(rjson)
library(data.table)
json <- fromJSON(file = "data/REACTOME_SUMOYLATION.v2022.1.Hs.json")
xgene <- json[[1]]$geneSymbols
save(xgene, file = "data/xgene.RData")

# 2.amigo
xgene <- fread("data/xgene.csv", sep = ":", header = F) %>% pull(V2)
suppressPackageStartupMessages(library(clusterProfiler))
conflict_prefer("first", "dplyr",quiet = T)
xgene <- bitr(geneID = xgene, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
xgene %<>% pull(2)
length(xgene)

# 3. txt file
xgene <- fread("data/Circadian Clock pathcards.txt",header = F) %>% pull(1)

# ---------------------------xgene与表达谱取交集------------------------------
load("data/train_data.RData")
length(xgene)
setdiff(xgene,rownames(train_data$data_exprs))
xgene <- intersect(rownames(train_data$data_exprs), xgene)
length(xgene)
save(xgene, file = str_glue("{od}/xgene.RData"))
