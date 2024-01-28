od <- "results/3.WGCNA"
source("src/0.config.R")
suppressWarnings(dir.create(od,recursive=TRUE))
load("data/train_data.RData")
load('results/2.DEG/com_gene.RData')
exprs_matrix <- train_data$tumor_exprs
head(exprs_matrix)
range(exprs_matrix)
# 做富集分析，后续将其作为一种临床信息的一种表型用于分析
library(GSVA)
ssgsea_res <- GSVA::gsva(
    expr = exprs_matrix %>% as.matrix(),
    gset.idx.list = list(com_gene),
    kcdf="Gaussian",
    method = "ssgsea",
    ssgsea.norm = TRUE,
    verbose = T,
    parallel.sz = 50
)
df=as.data.frame(ssgsea_res) %>% t() %>% as.data.table(keep.rownames = "Sample") %>% rename(Score=V1)
fwrite(df,paste0(od, "/ssgsea_score.csv"))

# 运行WGCNA，获得xgene相关基因
pheno <- as.data.frame(df)  %>% rename(sample=Sample)
head(pheno)

# install.packages('WGCNA')
# BiocManager::install('impute')
library(WGCNA)
library(conflicted)
conflicted::conflict_prefer("filter","dplyr")
qc_res <- wgcna_qc(
    exp = exprs_matrix,
    pheno = pheno,
    method = "average", 
	# cutHeight = 120,
	cutMad = 0.01, 
    od = paste0(od, "/QC"), width = 9, height = 6
)
# Before qc: nSample = 363; nGenes = 59427
# top 25% mad is 0.48
# top 50% mad is 0.00
# top 75% mad is 0.00; input cutMad is 0.01; max_mad is 0.01
#  Flagging genes and samples with too many missing values...
#   ..step 1
# 没有缺失值，不进行删除
# After qc: nSample = 363; nGenes = 26875
softpower_res <- wgcna_picksoftpower(exp = qc_res$use_exp, od = paste0(od, "/res"), net_type = "unsigned", Rsquared_cut = 0.8, width = 8, height = 6)
# The pick power: 14
# 根据选定参数构建WGCNA网络，生成ClusterDendrogram.pdf和ModuleTrait.pdf
wgcna_net <- wgcna_built_net(exp = qc_res$use_exp, pheno = qc_res$use_pheno, power = softpower_res$power, od = paste0(od, "/res"), height = 8, width = 6)
WGCNA::enableWGCNAThreads(nThreads = 8)
WGCNA::WGCNAnThreads()

write.table(wgcna_net$merged_infor, file = str_glue("{od}/res/merged_infor.txt"), quote = FALSE, sep = "\t", col.names = TRUE, row.names = F)

# 根据选定参数筛选WGCNA网络模块基因（核心基因），生成Figure_{select_pheno}_{select_module}_cor.pdf
merged_infor <- fread(str_glue("{od}/res/merged_infor.txt"),data.table=F)
colnames(merged_infor)
library(WGCNA)
library(conflicted)
conflict_prefer("filter", "dplyr")
head(merged_infor)
modulegene_res <- wgcna_select_modulegene(
    merged_infor = merged_infor, pheno_module_list = list(module = c("brown",'green'), 
    pheno =rep("Score",2)),
    od = paste0(od, "/colormdl") , MM = 0.6, GS = 0.3
)
# 在GS为0.3并且MM为0.6时，根据性状Score和模块brown，共筛选得到218个模块基因
# 在GS为0.3并且MM为0.6时，根据性状Score和模块green，共筛选得到34个模块基因

hub_gene <- flatten_chr(modulegene_res) %>% unique()
length(hub_gene)#252
fwrite(data.table(hub_gene=hub_gene),paste0(od,"/hub_gene.csv"))

DEG <- fread("output/SC02.score/CR_Score_Group_DEG_sig.csv") %>% pull(1)
WGCNA <- fread("results/3.WGCNA/hub_gene.csv") %>% pull(1)

library(ggvenn)
ggvenn(list(DEGs=DEG, WGCNA=WGCNA),show_percentage=F)
ggsave(paste0(od,'/venn.pdf'),width = 6,height = 6)

com_gene=intersect(DEG, WGCNA)
fwrite(data.table(com_gene=com_gene),paste0(od,"/com_gene2.csv"))
