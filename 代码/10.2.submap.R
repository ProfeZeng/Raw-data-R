conflict_prefer("filter", "dplyr")
# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")
od <- "results/10.treatment/"
dir.create(od)

valid_exp <- train_data$tumor_exprs
str(train_model_res,1)
train_model_res$Group %>% head()

# risk score group
com_sample = intersect(colnames(valid_exp),rownames(train_model_res$Group))
cls <- train_model_res$Group[com_sample,]


valid_exp=valid_exp[,com_sample]
rownames(valid_exp) %<>% str_replace_all('-','_')
colnames(valid_exp) %<>% str_replace_all('-','_')

valid_exp <- as.data.table(valid_exp,keep.rownames = 'NAME')
valid_exp %<>% mutate(description = '',.after='NAME')
fwrite(valid_exp,paste0(od,'/A.gct'),sep='\t')
nrow(valid_exp) # 59427
length(cls) # 363
fwrite(data.table(cls=cls),paste0(od,'/A.cls'),sep='\t')
#  --------------- submap ---------------
source('src/submap.R')
submap.main(
  input.data.A=paste0(od,'/A.gct'),
  input.data.B=paste0(od,'/skcm.immunotherapy.for.SubMap.gct'),
  input.cls.A=paste0(od,'/A.cls'),
  input.cls.B=paste0(od,'/skcm.immunotherapy.for.SubMap.cls'),
  output.filename="SubMap",
  ntag=100,
  nperm=50,
  nperm.fisher=1000,
  weighted.score.type=1,
  null.dist="pool",
  p.corr="FDR",
  clust.row=1,
  clust.col=1,
  nom.p.mat="F",
  create.legend="T",
  rnd.seed=1314
)
