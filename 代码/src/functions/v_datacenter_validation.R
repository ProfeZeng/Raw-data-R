#' @TODO 根据指定的癌型，循环遍历数据进行验证
#' @title 循环遍历数据进行验证(PCA方法以及普通的方法)
#' @description
#' @param signature_coef 数据框，包含两列signature和coef，模型的signature以及对应的系数coef；如果是PCA方法，则提供的是用于PCA分析的基因，基因的coef可以随意设置，不影响结果
#' @param od 结果输出路径
#' @param cancer 字符串，TCGA癌症简写，比如BRCA，一次只能填写一个
#' @param outcome_list 向量，对哪包含的临床随访信息进行遍历
#' @param edition 字符串，版本，默认v1
#' @param time ROC曲线生存时间节点，单位为年
#' @param best_cut 是否最优分组，默认FALSE
#' @param prognostic_independence 逻辑向量，是否进行预后独立性分析，默认FALSE
#' @param factors 字符串向量，用于独立性分析的临床特征,默认为c("Stage", "Age", "Gender")，不需要单独写Score
#' @param method 字符向量，默认normal，如果使用PCA的方式，则这里写成"PCA"
#' @param cohort_list 字符向量，验证集队列向量，默认NULL；如果指定的话，则只对指定队列进行验证
#' @param center 逻辑值，PCA分析的参数，默认为F（异常点多的情况下建议用F）
#' @return *list*
#' @examples
#' @author *CY*
#'
v_datacenter_validation <- function(signature_coef = NULL, od = getwd(),
                                    cancer = NULL, cox_tibble = NULL,no_roc_pheatmap = FALSE,
                                    outcome_list = c("OS", "PFS", "DFS", "DFI", "DSS", "PFI", "MFS", "RFS", "BCR"),
                                    edition = "v1", time = c(1, 3, 5), best_cut = FALSE,
                                    prognostic_independence = FALSE, factors = c("Stage", "Age", "Gender"),
                                    method = "normal", cohort_list = NULL, center = FALSE) {
    # 检查目录
    mkdir(od)
    if (is.null(cohort_list)) {
        valid_filelist <- readxl::read_xlsx(str_glue("{DATA_DIR}/GEO/idata_{edition}.xlsx")) %>%
            dplyr::filter(`癌症名称` == cancer) %>%
            dplyr::mutate(filepath = str_glue("{DATA_DIR}/{`来源数据库`}/{edition}/{`数据集名称`}/")) %>%
            dplyr::pull(filepath)
    } else {
        valid_filelist <- readxl::read_xlsx(str_glue("{DATA_DIR}/GEO/idata_{edition}.xlsx")) %>%
            dplyr::filter(`数据集名称` %in% cohort_list) %>%
            dplyr::mutate(filepath = str_glue("{DATA_DIR}/{`来源数据库`}/{edition}/{`数据集名称`}/")) %>%
            dplyr::pull(filepath)
    }
    message("exsit valid num is ", length(valid_filelist))
    # 模型遍历
    valid_filepath=valid_filelist
    valid_model_res <- map_dfr(valid_filelist, function(valid_filepath) {
        dataset <- str_split(string = valid_filepath, pattern = "/",simplify = T)[7]
        load(str_glue("{valid_filepath}/clinical.RData"))
        load(str_glue("{valid_filepath}/expression.RData"))
        # 提取表达谱
        data_exprs <- expression
        # 别名替换
        previous_name <- previous_name %>% dplyr::filter(PreviousORAliasSymbols %in% signature_coef$signature)
        if (nrow(previous_name) > 0) {
            for (i in 1:nrow(previous_name)) {
                if (any(rownames(data_exprs) == previous_name[i, "PreviousORAliasSymbols"]) & (sum(rownames(data_exprs) == previous_name[i, "ApprovedSymbol"]) < 1)) {
                    rownames(data_exprs)[rownames(data_exprs) == previous_name[i, "PreviousORAliasSymbols"]] <- previous_name[i, "ApprovedSymbol"]
                }
            }
        }
        data_clinical <- clinical %>%
            dplyr::filter(Sample %in% colnames(data_exprs), Sample.Type=="Tumor") %>%
            dplyr::rename(sample = Sample)
        if (any(colnames(data_clinical) %in% "Age")) {
            data_clinical <- clinical_pro_v2(input_clin = data_clinical, agecol = "Age")
        }
        data_exprs[is.na(data_exprs)] <- 0
        common_sample <- intersect(colnames(data_exprs), data_clinical$sample)
        common_gene <- intersect(rownames(data_exprs), signature_coef$signature)
        data_exprs <- data_exprs[, common_sample]
        next_step <- 0
        if (!str_detect(method, "PCA") & length(common_gene) >= 1) {
            exprs <- data_exprs
            next_step <- 1
        } else if (method == "PCA1" & length(common_gene) >= 5) {
            pca_res <- prcomp(x = data_exprs[common_gene, common_sample], retx = F, scale = T, center = center)
            exprs <- pca_res$rotation %>%
                as.data.frame() %>%
                mutate(PC1PC2 = PC1 + PC2) %>%
                dplyr::select(PC1PC2) %>%
                t()
            signature_coef <- data.frame(signature = "PC1PC2", coef = 1)
            next_step <- 1
        } else if (method == "PCA2" & length(common_gene) >= 5) {
            cox_tibble <- cox_tibble[common_gene, ]
            HR_high_genes <- cox_tibble %>%
                dplyr::filter(HR > 1) %>%
                rownames()
            HR_low_genes <- cox_tibble %>%
                dplyr::filter(HR < 1) %>%
                rownames()
            p_load(FactoMineR)
            pcaH <- PCA(X=data_exprs[HR_high_genes, ] %>% t(),graph = FALSE, ncp = 2)
            pcaL <- PCA(X=data_exprs[HR_low_genes, ] %>% t(),graph = FALSE, ncp = 2)
            score <- pcaH$ind$contrib[, 1] + pcaH$ind$contrib[, 2] - pcaL$ind$contrib[, 1] - pcaL$ind$contrib[, 2]
            exprs <- data.frame(PC1PC2 = score) %>% t()
            signature_coef <- tibble(signature = "PC1PC2", coef = 1)
            next_step <- 1
        }
        if (next_step == 1) {
            model_res <- map_dfr(outcome_list, function(survival_outcome) {
                if (any(colnames(data_clinical) %in% str_glue("{survival_outcome}.Time")) & any(colnames(data_clinical) %in% str_glue("{survival_outcome}.Status"))) {
                    source("src/functions/v_modelscore_km_roc.R",local = TRUE)
                    res <- v_modelscore_km_roc(
                        signature_coef = signature_coef,
                        od = str_glue("{od}/{survival_outcome}/"),
                        exp = exprs,
                        clin = data_clinical,
                        dataset = dataset, 
                        time = time,
                        no_roc_pheatmap = no_roc_pheatmap,
                        best_cut = best_cut,
                        timecol = str_glue("{survival_outcome}.Time"), 
                        statuscol = str_glue("{survival_outcome}.Status"),
                        xlab = str_c(survival_outcome, " days")
                    )
                    if (!is.null(res)) {
                        return_res <- res$infor %>% mutate(survival_outcome = survival_outcome, cancer = cancer, common_gene_num = length(common_gene))
                        if (prognostic_independence) {
                            source("src/functions/v_uni_multi_cox.R", local = TRUE)
                            uni_multi_cox(
                                od = str_glue("{od}/{survival_outcome}/"), infor = res$Group %>%
                                    rownames_to_column("sample") %>% dplyr::rename(Score = Group) %>%
                                    merge(., data_clinical), factors = c(factors, "Score"), dataset = dataset, coxp = 1, w = 10,
                                timecol = str_glue("{survival_outcome}.Time"), statuscol = str_glue("{survival_outcome}.Status")
                            )
                        }
                    } else {
                        return_res <- data.frame()
                    }
                } else {
                    return_res <- data.frame()
                }
                return(return_res)
            })
            return(model_res)
        } else {
            return(data.frame())
        }
    })
    return(valid_model_res)
}

# source("/Pub/Users/cuiye/RCodes/UserCode/newlover/plotout.R")
# source("/Pub/Users/cuiye/RCodes/UserCode/newlover/mkdir.R")
# # previous_name <- read.delim(file = "/Pub/Users/cuiye/database/HUGO_Approved_Alias_symbols.txt", header = TRUE, check.names = FALSE)
# # save(previous_name, file = "/Pub/Users/cuiye/database/previous_name.RData")
# load("/Pub/Users/cuiye/database/previous_name.RData")
