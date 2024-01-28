# -------------------导入函数-------------------
source("src/0.config.R")
using(tidyverse)
# -------------------设置结果路径-------------------
od <- "data"
using(arrow)
df=arrow::read_ipc_file("results/1.prepare/TCGA_LIHC_mRNA.arrow") %>% column_to_rownames('gene_name') 
df=log2(df+1)
df[1:5,1:5]
data_pdata = data.frame(Sample=colnames(df),Group=ifelse(str_detect(colnames(df),'^.{13}01A'),'Tumor','Normal'))
train_data=list(data_exprs=df,data_pdata=data_pdata,tumor_exprs=dplyr::select(df,matches('^.{13}01A')))
train_data$tumor_exprs[1:5,1:5]
dim(train_data$tumor_exprs)
str(train_data,1)
# -------------------其它临床指标-------------------
data_clinical <- readxl::read_xlsx(paste0(Data_Center,'/TCGA_secondary/TCGA-CDR-SupplementalTableS1.xlsx'),sheet = 'TCGA-CDR') %>% 
    dplyr::filter(type %in% parameter_list$cancer[1]) %>%
    dplyr::mutate(OS.Time = !!sym(str_c(parameter_list$train_survival_outcome, ".time")), 
        OS.Status = !!sym(parameter_list$train_survival_outcome)) %>%
    dplyr::mutate(OS.Time = as.numeric(OS.Time), OS.Status = as.numeric(OS.Status)) %>%
    dplyr::filter(OS.Time > 0 & OS.Status != "" & OS.Status != "NA") %>% 
    dplyr::mutate(Age=age_at_initial_pathologic_diagnosis,Sex=gender,
        ,Stage=case_when(str_detect(ajcc_pathologic_tumor_stage,'Stage I[A-Ea-e]$')~'I', #/II
        str_detect(ajcc_pathologic_tumor_stage,'Stage II[A-Ea-e]$')~'II',#
        str_detect(ajcc_pathologic_tumor_stage,'Stage III[A-Ea-e]$')~'III',#
        str_detect(ajcc_pathologic_tumor_stage,'Stage IV[A-Ea-e]$')~'IV') #III/IV
        ,Sample=paste0(bcr_patient_barcode,'-01A')
        ,Age.Detail = Age
        ,sample=Sample
        ,Age=ifelse(Age.Detail>=60,"Age>=60","Age<60")
        )

data_clinical$Race <- ifelse(str_detect(data_clinical$race,"\\["),NA,data_clinical$race)
data_clinical$Grade <- ifelse(str_detect(data_clinical$histological_grade,"\\["),NA,data_clinical$histological_grade)
data_clinical$tumor_status <- ifelse(str_detect(data_clinical$tumor_status,"\\["),NA,data_clinical$tumor_status)

train_data$data_clinical <- data_clinical
# -------------------保留临床随访信息齐全和表达量均有的肿瘤样本-------------------
keep_tumor_sample <- intersect(train_data$data_clinical$Sample, colnames(train_data$tumor_exprs))
train_data$data_clinical <- train_data$data_clinical %>% dplyr::filter(Sample %in% keep_tumor_sample)
train_data$tumor_exprs <- train_data$tumor_exprs[, keep_tumor_sample]
save(train_data, file = str_glue("{od}/train_data.RData"))

write.table(train_data$tumor_exprs, file = str_glue("{od}/SupplementaryTable_train_tumor_exprs.txt"), quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
write.table(train_data$data_clinical, file = str_glue("{od}/SupplementaryTable_train_data_clinical.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# -------------------统计样本的临床特征-------------------
stat <- train_data$data_clinical %>%
    dplyr::select(Sample, all_of(parameter_list$clinical_feature_list)) %>%
    pivot_longer(
        cols = -Sample,
        names_to = "type",
        values_to = "value"
    ) %>%
    group_by(type, value) %>%
    summarise(sample_num = n())
write.table(stat, file = str_glue("{od}/SupplementaryTable_train_data_clinical_stat.txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

load("data/train_data.RData")
fwrite(as.data.table(train_data$tumor_exprs,keep.rownames = "Sample"),"data/TCGA_Tumor.csv.gz")
load("results/5.model/signature/train_model_res.RData")
fwrite(as.data.table(train_model_res$Group ,keep.rownames = "Sample"),"data/TCGA_Group.csv")
