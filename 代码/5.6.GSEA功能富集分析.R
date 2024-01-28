rm(list=ls())
source("")
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)

# 设置结果路径
od <- "results/5.model/GSEA"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("results/5.model/DEG/deg_gene.RData")

# 排序
geneList <- deg_res$nrDEG$logFC
names(geneList) <- rownames(deg_res$nrDEG)
geneList <- sort(geneList, decreasing = TRUE)
# 富集分析
library(clusterProfiler)
library(enrichplot)
Data_Center
geneSets <- read.gmt(str_glue("{Data_Center}/Gene_Sets/MsigDB/msigdb.v2022.1.Hs.symbols.gmt"))

lapply(c("KEGG", "GOBP","HALLMARK"), function(database) {
    gsea_result <- GSEA(geneList,
        TERM2GENE = geneSets %>% filter(grepl(pattern = database, x = `term`, ignore.case = T)),
        eps = 0, pvalueCutoff = 0.05, pAdjustMethod = "BH"
    )
    if (nrow(as.data.frame(gsea_result)) > 0) {
        write.table(gsea_result, file = str_glue("{od}/SupplementaryTable_{database}_term.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
        # 画图
        if (nrow(as.data.frame(gsea_result)) >= 20) {
            n <- 20
        } else {
            n <- nrow(as.data.frame(gsea_result))
        }
        lapply(1:n, function(x) {
            p <- gseaplot2(gsea_result, x, color = "red", pvalue_table = T)
            plotout(p = p, h = 6, w = 10, od = str_glue("{od}/{database}/"), name = rownames(as.data.frame(gsea_result))[x])
            return(NULL)
        })
    }
    return(NULL)
})

# 挑选结果画图
if (FALSE) {
    library(clusterProfiler)
    library(enrichplot)
    database <- "KEGG"
    geneSets <- read.gmt("/Pub/Data/Data_Center/GeneSet/MsigDB/v7.4/msigdb.v7.4.symbols.gmt")
    gsea_result <- GSEA(geneList,
        TERM2GENE = geneSets %>% filter(grepl(pattern = database, x = `term`, ignore.case = T)),
        eps = 0, pvalueCutoff = 0.05, pAdjustMethod = "BH"
    )
    path <- c(
        "KEGG_ECM_RECEPTOR_INTERACTION", "KEGG_CELL_ADHESION_MOLECULES_CAMS", "KEGG_CALCIUM_SIGNALING_PATHWAY",
        "KEGG_REGULATION_OF_ACTIN_CYTOSKELETON", "KEGG_MAPK_SIGNALING_PATHWAY"
    )
    p <- gseaplot2(
        x = gsea_result, geneSetID = path,
        pvalue_table = TRUE,
        color = color_fun1[1:5],
        ES_geom = "line"
    ) + xlab("High vs Low")
    plotout(od = od, name = "select_KEGG", w = 12, h = 8, p = p, plot_tif = FALSE)
}
