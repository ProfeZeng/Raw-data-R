# ==============================================================================
# 1.安装和加载R包
# ==============================================================================
using(data.table, tidyverse, CellChat, NMF, ggalluvial, patchwork, arrow, reticulate, ComplexHeatmap)
reticulate::use_python("/opt/homebrew/Caskroom/mambaforge/base/envs/SC/bin/python")
output_dir <- "output/SC05.CellChat"
org <- "human"
id <- "High"
group_names <- c("CR_High", "CR_Low")

psave <- function(filename, w = 12, h = 6, plot = FALSE) {
    if (is.object(p)) {
        print(p)
    }
    plot <- recordPlot()
    pdf(file = filename, onefile = T, width = w, height = h)
    replayPlot(plot)
    dev.off()
}

# ==============================================================================
# 5. 合并对象
# ==============================================================================
cc1 <- readRDS(str_glue("{output_dir}/cc_{group_names[1]}.rds"))
cc2 <- readRDS(str_glue("{output_dir}/cc_{group_names[2]}.rds"))
cc_list <- list(cc1, cc2)
names(cc_list) <- group_names
cc <- mergeCellChat(cc_list, add.names = group_names)

# ==============================================================================
# 6.
# ==============================================================================
p1 <- compareInteractions(cc, show.legend = FALSE, group = c(1, 2))
p2 <- compareInteractions(cc, show.legend = FALSE, group = c(1, 2), measure = "weight")
p1 + p2
psave(str_glue("{output_dir}/Overview_number_strength.pdf"), w = 6, h = 4)

# ==============================================================================
# 7. 数量与强度差异网络图
# ==============================================================================
par(mfrow = c(1, 2))
netVisual_diffInteraction(cc, weight.scale = T)
netVisual_diffInteraction(cc, weight.scale = T, measure = "weight")
psave(str_glue("{output_dir}/Diff_number_strength_net.pdf"))
netVisual_heatmap(cc, measure = "count")
psave(str_glue("{output_dir}/Diff_number_strength_heatmap_count.pdf"))

# ==============================================================================
# 8. 数量与强度差异热图
# ==============================================================================
netVisual_heatmap(cc, measure = "weight")
psave(str_glue("{output_dir}/Diff_number_strength_weight.pdf"))
par(mfrow = c(1, 2))
weight.max <- getMaxWeight(cc_list, attribute = c("idents", "count"))
for (i in 1:length(cc_list)) {
    netVisual_circle(cc_list[[i]]@net$count,
        weight.scale = T, label.edge = F,
        edge.weight.max = weight.max[2], edge.width.max = 12,
        title.name = paste0("Number of interactions - ", names(cc_list)[i])
    )
}
psave(str_glue("{output_dir}/Counts_Compare_net.pdf"))

# ==============================================================================
# 指定细胞互作数量对比网络图
# ==============================================================================

par(mfrow = c(1, 2))
s.cell <- c("CD4+ T cells", "CD8+ T cells", "Monocytes")
count1 <- cc_list[[1]]@net$count[s.cell, s.cell]
count2 <- cc_list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1,
    weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12,
    title.name = paste0("Number of interactions-", names(cc_list)[1])
)
netVisual_circle(count2,
    weight.scale = T, label.edge = T, edge.weight.max = weight.max, edge.width.max = 12,
    title.name = paste0("Number of interactions-", names(cc_list)[2])
)
# save as Counts_Compare_select.pdf 10*6.5

# ==============================================================================
# 保守和特异性信号通路的识别与可视化 ## 通路信号强度对比分析
# ==============================================================================

p1 <- rankNet(cc, mode = "comparison", stacked = TRUE, do.stat = TRUE)
p2 <- rankNet(cc, mode = "comparison", stacked = FALSE, do.stat = TRUE)
p1 + p2
psave(str_glue("{output_dir}/Compare_pathway_strengh.pdf"), w = 10, h = 10)

# ==============================================================================
# 流行学习识别差异信号通路
# ==============================================================================

cc <- computeNetSimilarityPairwise(cc, type = "functional")
cc <- netEmbedding(cc, type = "functional") # Python
cc <- netClustering(cc, type = "functional", do.parallel = FALSE)

netVisual_embeddingPairwise(cc, type = "functional", label.size = 3.5)
psave(str_glue("{output_dir}/functional.pdf"), w = 4, h = 3)

netVisual_embeddingPairwiseZoomIn(cc, type = "functional", nCol = 3)
psave(str_glue("{output_dir}/functional_zoom.pdf"), w = 12, h = 10)

rankSimilarity(cc, type = "functional") + ggtitle("functional similarity of pathway")
psave(str_glue("{output_dir}/Pathway_functional_similarity.pdf"), w = 9, h = 6)

cc <- computeNetSimilarityPairwise(cc, type = "structural")
cc <- netEmbedding(cc, type = "structural")
cc <- netClustering(cc, type = "structural", do.parallel = FALSE)

netVisual_embeddingPairwise(cc, type = "structural", label.size = 3.5)
psave(str_glue("{output_dir}/structural.pdf"), w = 4, h = 3)

netVisual_embeddingPairwiseZoomIn(cc, type = "structural", nCol = 2)
psave(str_glue("{output_dir}/structural_zoom.pdf"), w = 12, h = 10)

rankSimilarity(cc, type = "structural") + ggtitle("Structural similarity of pathway")
psave(str_glue("{output_dir}/Pathway_structural_similarity.pdf"), w = 9, h = 6)

# ==============================================================================
# 3.5 细胞信号模式对比
# ==============================================================================

dir.create(str_glue("{output_dir}/Compare_Signal"), recursive = TRUE)
pathway_union <- union(cc_list[[1]]@netP$pathways, cc_list[[2]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(cc_list[[1]], title = group_names[1], pattern = "all", signaling = pathway_union, width = 12, height = 28)
ht2 <- netAnalysis_signalingRole_heatmap(cc_list[[2]], title = group_names[2], pattern = "all", signaling = pathway_union, width = 12, height = 28)
ht1 + ht2
psave(str_glue("{output_dir}/Compare_Signal/Signal_all.pdf"), w = 12, h = 13)

ht1 <- netAnalysis_signalingRole_heatmap(cc_list[[1]], title = group_names[1], pattern = "outgoing", signaling = pathway_union, width = 12, height = 28)
ht2 <- netAnalysis_signalingRole_heatmap(cc_list[[2]], title = group_names[2], pattern = "outgoing", signaling = pathway_union, width = 12, height = 28)
ht1 + ht2
psave(str_glue("{output_dir}/Compare_Signal/Signal_outgoing.pdf"), w = 12, h = 13)

ht1 <- netAnalysis_signalingRole_heatmap(cc_list[[1]], title = group_names[1], pattern = "incoming", signaling = pathway_union, width = 12, height = 28)
ht2 <- netAnalysis_signalingRole_heatmap(cc_list[[2]], title = group_names[2], pattern = "incoming", signaling = pathway_union, width = 12, height = 28)
ht1 + ht2
psave(str_glue("{output_dir}/Compare_Signal/Signal_incoming.pdf"), w = 12, h = 13)

# ==============================================================================
# 3.6 特定信号通路的对比
# ==============================================================================
setdiff(cc_list[[2]]@netP$pathways,cc_list[[1]]@netP$pathways)

com_pathways <- intersect(cc_list[[1]]@netP$pathways, cc_list[[2]]@netP$pathways)
dir.create(str_glue("{output_dir}/Compare_Pathways_Net"), recursive = TRUE)
for (x in com_pathways) {
    weight.max <- getMaxWeight(cc_list, slot.name = c("netP"), attribute = x)
    par(mfrow = c(1, 2), xpd = TRUE)

    for (y in 1:length(cc_list)) {
        netVisual_aggregate(cc_list[[y]],
            signaling = x,
            layout = "circle",
            edge.weight.max = weight.max[1],
            edge.width.max = 10,
            signaling.name = paste(x, group_names[y])
        )
    }
    psave(str_glue("{output_dir}/Compare_Pathways_Net/Compare_{x}_Net.pdf"), w = 12, h = 6)
}

dir.create(str_glue("{output_dir}/Compare_Pathways_Chord"), recursive = TRUE)
for (x in com_pathways) {
    weight.max <- getMaxWeight(cc_list, slot.name = c("netP"), attribute = x)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (y in 1:length(cc_list)) {
        netVisual_aggregate(cc_list[[y]],
            signaling = x,
            layout = "chord",
            pt.title = 3,
            title.space = 0.05,
            vertex.label.cex = 0.6,
            edge.weight.max = weight.max[1],
            edge.width.max = 10,
            signaling.name = paste(x, group_names[y])
        )
    }
    psave(str_glue("{output_dir}/Compare_Pathways_Chord/Compare_{x}_Chord.pdf"), w = 9, h = 6)
}

dir.create(str_glue("{output_dir}/Compare_Pathways_Heatmap"), recursive = TRUE)

for (x in setdiff(com_pathways, c("PECAM1", "OX40", "PVR", "TIGIT", "CALCR"))) {
    par(mfrow = c(1, 2), xpd = TRUE)
    ht <- list()
    for (y in 1:length(cc_list)) {
        # y=1
        ht[[y]] <- netVisual_heatmap(
            cc_list[[y]],
            signaling = x,
            title.name = paste(x, "signaling ", group_names[y])
        )
    }
    p <- purrr::reduce(ht, `+`)
    print(p)
    psave(str_glue("{output_dir}/Compare_Pathways_Heatmap/{x}.pdf"))
}

# ==============================================================================
# 气泡图展示High组上调或下调的配体受体对
# ==============================================================================

dir.create(str_glue("{output_dir}/Compare_LR_regulated"), recursive = TRUE)
for (x in setdiff(com_pathways, c("CXCL", "PECAM1", "LCK", "OX40", "PVR", "TIGIT", "CALCR", "GRN", "CDH5", "BMP", "MPZ")))
{
    p <- netVisual_bubble(cc,
        # sources.use = c(4, 5), targets.use = c(1, 2, 3, 6),
        comparison = c(1, 2), signaling = x,
        max.dataset = 1, title.name = str_glue("Increased signaling in {group_names[1]}"), angle.x = 45, remove.isolate = T
    )

    psave(str_glue("{output_dir}/Compare_LR_regulated/{x}_Increased.pdf"), plot = p, w = 18, h = 6)
    p <- netVisual_bubble(cc,
        # sources.use = c(4, 5), targets.use = c(1, 2, 3, 6),
        comparison = c(1, 2), signaling = x,
        min.dataset = 1, title.name = str_glue("Decreased signaling in {group_names[1]}"), angle.x = 45, remove.isolate = T
    )
    psave(str_glue("{output_dir}/Compare_LR_regulated/{x}_Decreased.pdf"), plot = p, w = 18, h = 6)
}

saveRDS(cc, str_glue("{output_dir}/cc.rds"))
cc <- readRDS(str_glue("{output_dir}/cc.rds"))
