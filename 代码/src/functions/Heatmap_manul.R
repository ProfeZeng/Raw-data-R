#' @TODO 组间分类绘制热图
#' @title 组间分类绘制热图
#' @description 利用`ComplexHeatmap`绘制热图，默认对行进行标准化，标准化范围 默认c(-1,1).
#' @param data_input 输入数据，数据框，比如表达谱
#' @param color_used 组间分类使用配色，如果为`NULL`则默认Set1配色
#' @param group_infor 分组信息，可以是单列或者多列。
#' 格式要求：样本在行名，用第一列作为配色列或者说分组依据，例如聚类分组。
#' group可以是任意分组信息
#' @param Colored_OtherInfor *logical value* or *character string* ，是否对分组信息中的其他信息进行染色，默认使用Paired调色板，
#' 如果是颜色字符串，则使用自定义颜色对其他分组信息按顺序染色，
#' @param saveplot *logical value*，是否保存图片
#' @param output_dir *character string*，文件输出路径，最好以/结尾
#' @param var_name *character string*，成图文件命名字段，字符串向量
#' @param width_used *numric*，图宽
#' @param height_used *numric*，图高
#' @param cluster_name *character string*，分组名字可以是’cluster'等字段,已弃用，默认使用group_infor第一列做为分组依据
#' @param heatmap_name *character string*，热图colorbar name，默认为NULL
#' @param DoWilcox.test *logical value*，是否对热图主体matirx 在组间使用秩和检验，检验显著性差异
#' @param show_rownames *logical value*，是否显示行名
#' @return  ggplot2对象，热图
#' @usage
#' a <- Heatmap_manul(
#'    data_input = expr[geneList, ],
#'    color_used = RColorBrewer::brewer.pal(3, "Set2"),
#'    cluster_name = 'cluster',
#'    group_infor = group_data %>% rename(group = 1),
#'    saveplot = T,
#'    output_dir = output_dir,
#'    var_name = "DLBC",
#'    width_used = 10,
#'    height_used = 12
#' )
#' @usage 表达谱示例
#' > expr[geneList, ][1:3,1:3]
#'          GSM775979 GSM775980  GSM775981
#' ATOX1    0.4598178 0.4990092 0.22533183
#' C1orf122 0.2583019 0.1842661 0.26094508
#' CENPE    0.2950682 0.3843891 0.06562005
#' @usage 分组信息示例
#' > group_data %>% rename(group = 1) %>% .[1:3,1:3]
#'           group     Age gender
#' GSM775979     A  Age<60      M
#' GSM775980     A Age>=60      F
#' GSM775981     B  Age<60      F
#' @export
#' @author *WYK*
Heatmap_manul <- function(data_input = NULL, color_used = NULL, group_infor = NULL, Colored_OtherInfor = F,
    saveplot = T, output_dir = "./", var_name = NULL, width_used = 9, height_used = 9, heatmap_name = " ",
    cluster_name = NULL, show_rownames = T, DoWilcox.test = F, rownames_fontsize = 7) {
    conflicts_prefer(dplyr::select)
    if (is.null(var_name)) {
        var_name <- paste0("_", paste0(sample(letters, 4), collapse = "", sep = ""))
    }

    if (!is.null(cluster_name)) {
        warning(stringr::str_glue('在function {crayon::blue("Heatmap_manul")}中，参数{crayon::bold("cluster_name")} 已弃用,并且使用{crayon::bold("group_infor")}中第一列做分组依据.'))
    }
    suppressPackageStartupMessages(library(ComplexHeatmap))
    suppressPackageStartupMessages(library(tidyverse))
    suppressPackageStartupMessages(library(cowplot))
    suppressPackageStartupMessages(library(ggpubr))
    cluster_name <- colnames(group_infor)[1]

    sample_common <- intersect(rownames(group_infor), colnames(data_input))
    group_infor <- group_infor[sample_common, , drop = F] %>% dplyr::rename("group" = 1)

    data_input <- data_input[, sample_common]

    # if (ncol(group_infor) >= 2) {
    #     clinical_names <- base::setdiff(colnames(group_infor), "group")
    #     colnames_chr <- map_chr(clinical_names, function(x) {
    #         print(x)
    #         fisher_test_res <- group_infor %>%
    #             dplyr::select(group, any_of(x)) %>%
    #             table() %>%
    #             fisher.test(x=.,simulate.p.value = TRUE, B = 1e7)
    #         fisher_p <- fisher_test_res$p.value
    #         pval_label <- cut(
    #             x = fisher_p, breaks = c(1, .05, .01, .001, .0001, 0),
    #             labels = c("****", "***", "**", "*", "")
    #         ) %>% as.character()

    #         x_pval_label <- paste0(x, " ", pval_label)
    #     })

    #     colnames(group_infor)[c(1:ncol(group_infor))[-which(colnames(group_infor) %in% "group")]] <- colnames_chr
    # }

    group_chara <- sort(unique(as.character(group_infor$group)))

    if (is.null(color_used)) {
        col_name <- RColorBrewer::brewer.pal(8, "Set1")[1:length(group_chara)]
        col_name <- col_name[seq_along(unique(group_infor %>% pull(group)))]
    } else {
        col_name <- color_used[seq_along(unique(group_infor %>% pull(group)))]
    }

    names(col_name) <- group_chara

    if (ncol(group_infor) >= 2) {
        col_anno <- HeatmapAnnotation(
            df = group_infor, # %>% .[, "group", drop = F]
            col = list(
                group = col_name
            ),
            annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
            annotation_name_side = "right",
            annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
        )
    } else {
        col_anno <- HeatmapAnnotation(
            df = group_infor, # %>% .[, "group", drop = F]
            col = list(
                group = col_name
            ),
            annotation_name_side = "right",
            annotation_label = cluster_name
        )
    }

    if (isTRUE(Colored_OtherInfor)) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")

        color_defined <- rep(RColorBrewer::brewer.pal(9, "Paired"), 10)
        # color_defined <- rep(RColorBrewer::brewer.pal(12,'Set3'),10)
        i <- 1

        color_in_heatmap <- lapply(clinical_names, function(x) {
            # x <- clinical_names[1]
            color_name <- group_infor %>%
                pull(x) %>%
                unique() %>%
                sort()

            color_used <- color_defined[i:(i + length(color_name) - 1)]
            names(color_used) <- color_name

            color_defined <- color_defined[-c(i:(i + length(color_name) - 1))]
            i <<- i + length(color_name)
            return(color_used)
        })

        names(color_in_heatmap) <- clinical_names

        color_in_heatmap$group <- col_name

        if (ncol(group_infor) >= 2) {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = color_in_heatmap,
                annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
                annotation_name_side = "right",
                annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
            )
        } else {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = list(
                    group = col_name
                ),
                annotation_name_side = "right",
                annotation_label = cluster_name
            )
        }
    }

    if (is.null(Colored_OtherInfor)) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")

        color_defined <- rep(RColorBrewer::brewer.pal(9, "Paired"), 10)
        # color_defined <- rep(RColorBrewer::brewer.pal(12,'Set3'),10)
        i <- 1

        color_in_heatmap <- lapply(clinical_names, function(x) {
            # x <- clinical_names[1]
            color_name <- group_infor %>%
                pull(x) %>%
                unique() %>%
                sort()

            color_used <- color_defined[i:(i + length(color_name) - 1)]
            names(color_used) <- color_name

            color_defined <- color_defined[-c(i:(i + length(color_name) - 1))]
            i <<- i + length(color_name)
            return(color_used)
        })

        names(color_in_heatmap) <- clinical_names
        color_in_heatmap$group <- col_name

        if (ncol(group_infor) >= 2) {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = color_in_heatmap,
                annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
                annotation_name_side = "right",
                annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
            )
        } else {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = list(
                    group = col_name
                ),
                annotation_name_side = "right",
                annotation_label = cluster_name
            )
        }
    }

    if (is.character(Colored_OtherInfor)) {
        clinical_names <- base::setdiff(colnames(group_infor), "group")

        color_defined <- rep(Colored_OtherInfor, 30)
        # color_defined <- rep(RColorBrewer::brewer.pal(12,'Set3'),10)
        i <- 1

        color_in_heatmap <- lapply(clinical_names, function(x) {
            # x <- clinical_names[1]
            color_name <- group_infor %>%
                pull(x) %>%
                unique() %>%
                sort()

            color_used <- color_defined[i:(i + length(color_name) - 1)]
            names(color_used) <- color_name

            color_defined <- color_defined[-c(i:(i + length(color_name) - 1))]
            i <<- i + length(color_name)
            return(color_used)
        })

        names(color_in_heatmap) <- clinical_names

        color_in_heatmap$group <- col_name

        if (ncol(group_infor) >= 2) {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = color_in_heatmap,
                annotation_legend_param = list(group = list(title = c(cluster_name, base::setdiff(colnames(group_infor), "group")))),
                annotation_name_side = "right",
                annotation_label = c(cluster_name, base::setdiff(colnames(group_infor), "group"))
            )
        } else {
            col_anno <- HeatmapAnnotation(
                df = group_infor, # %>% .[, "group", drop = F]
                col = list(
                    group = col_name
                ),
                annotation_name_side = "right",
                annotation_label = cluster_name
            )
        }
    }

    data_input_scaled <- t(apply(data_input, 1, function(x) scale(x, center = T, scale = T)))
    colnames(data_input_scaled) <- colnames(data_input)

    # data_input_scaled <- map_df(1:nrow(data_input), function(i) {
    #     z_scores <- (data_input[i, ] - mean(as.numeric(data_input[i, ]))) / sd(as.numeric(data_input[i, ]))
    #     return(z_scores)
    # })

    all(colnames(data_input_scaled) == group_infor %>% rownames())
    col_zscore <- circlize::colorRamp2(c(-2, 0, 2), c("#2266AC", "white", "#B2182E")) # c("#0a5aa5", "white", "firebrick3")

    # col_zscore <- circlize::colorRamp2(c(-2, 0, 2), c("#00441B", "white", "#40004B"))

    cluster_res <- Heatmap(
        matrix = data_input_scaled,
        name = heatmap_name,
        col = col_zscore,
        show_column_dend = F,
        show_column_names = F,
        show_row_names = show_rownames,
        cluster_rows = F,
        cluster_columns = F,
        show_row_dend = F,
        clustering_method_columns = "complete",
        # rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = rownames_fontsize),
        column_title = " ",
        top_annotation = col_anno,
        column_split = group_infor %>% .[, "group", drop = F],
        border = F
        # left_annotation = row_anno,
        # row_names_side = "right"
    )

    if (isTRUE(DoWilcox.test)) {
        group_wilcoxtest_p <- map_dbl(1:nrow(data_input), function(i) {
            # i <- 1
            if (length(unique(group_infor[, 1])) == 2) {
                data_used <- data.frame(value = data_input[i, ] %>% as.numeric(), group = group_infor[match(colnames(data_input), rownames(group_infor)), "group"])
                wilcox.test_res <- wilcox.test(value ~ group, data = data_used)
                pval <- wilcox.test_res$p.value

                return(pval)
            } else {
                data_used <- data.frame(value = data_input[i, ] %>% as.numeric(), group = group_infor[match(colnames(data_input), rownames(group_infor)), "group"])
                kruskal.test_res <- kruskal.test(value ~ group, data = data_used)
                pval <- kruskal.test_res$p.value

                return(pval)
            }
        })

        names(group_wilcoxtest_p) <- rownames(data_input)

        row_anno <- rowAnnotation(
            P.val = anno_text(case_when(
                between(group_wilcoxtest_p, 0.01, 0.05) ~ "*",
                between(group_wilcoxtest_p, 0.001, 0.01) ~ "**",
                between(group_wilcoxtest_p, 0.0001, 0.001) ~ "***",
                group_wilcoxtest_p < 0.0001 ~ "****",
                group_wilcoxtest_p > 0.05 ~ " "
            ),
            gp = gpar(fontsize = 8),
            location = 1,
            just = "right"
            )
        )

        cluster_res <- Heatmap(
            matrix = data_input_scaled,
            name = heatmap_name,
            col = col_zscore,
            show_column_dend = F,
            show_column_names = F,
            show_row_names = show_rownames,
            cluster_rows = F,
            cluster_columns = F,
            show_row_dend = F,
            clustering_method_columns = "complete",
            # rect_gp = gpar(col = "white", lwd = 1),
            row_names_gp = gpar(fontsize = rownames_fontsize),
            column_title = " ",
            top_annotation = col_anno,
            column_split = group_infor %>% .[, "group", drop = F],
            border = F,
            left_annotation = row_anno
            # row_names_side = "right"
        )
    }

    if (saveplot) {
        dir_now <- str_glue("{output_dir}")

        if (!dir.exists(dir_now)) {
            dir.create(dir_now, recursive = T)
        } else {
            message(str_c(dir_now, " is ready."))
        }

        pdf(
            file = "Figure_Heatmap.pdf",
            width = width_used, height = height_used
        )
        draw(cluster_res, merge_legend = T)
        dev.off()
    }
    return(cluster_res)
}