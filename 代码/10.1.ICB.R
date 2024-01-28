#---------------------------TIDE analysis------------------------------
rm(list=ls())
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)
# 设置结果路径
od <- "results/10.treatment"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
load("data/xgene.RData")
load("results/5.model/signature/train_model_res.RData")
tumor_exprs <- train_data$tumor_exprs

#ouput for TIDE
Expr <- t(apply(tumor_exprs, 1, function(x) x-mean(x)))
idmaps = clusterProfiler::bitr(geneID=rownames(Expr),fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = TRUE)
Expr2 = as_tibble(Expr,rownames = "SYMBOL") %>% 
    inner_join(x=idmaps,y=.,by="SYMBOL") %>% 
    dplyr::select(-SYMBOL)

fwrite(Expr2,file=str_glue("{od}/expression_4_TIDE.txt"),sep="\t")

# sehll
tidepy expression_4_TIDE.txt -o tide_res.tsv -c Other

# read tide res
tide_res = fread(str_glue('{od}/tide_res.tsv'))
tide_res[1:3,1:3]
plot_data = merge(tide_res,as_tibble(train_model_res$Group,rownames = "sample"),by.x='V1',by.y="sample")
plot_data %<>% dplyr::mutate(Responder=factor(Responder,levels = c("FALSE","TRUE")))

#boxplot
library(ggpubr)
plot_data$Responder
ggboxplot(plot_data,palette=color_fun1,order=c('TRUE','FALSE'),
    x = "Responder", y = "TIDE", fill = "Responder",
    xlab = "Responder", ylab = "TIDE Value",
    title = "", add = ""
) +
    stat_compare_means(label.y.npc = 1, label.x.npc = 0.4,size=5) +
    theme(plot.title=element_text(size=20,hjust=0.5),
        panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = "", 
        axis.title.x = element_text(size=14),
        axis.title.y = element_text( size=14),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(color = "black")
    )
ggsave(filename = str_glue("{od}/Figure_TID_boxplot.pdf"),width=5,height = 5)

ggboxplot(plot_data,palette=color_fun1,
    x = "Group", y = "TIDE", fill = "Group",
    xlab = "", ylab = "TIDE Value",
    title = "", add = ""
) +
    stat_compare_means(label.y.npc = 1, label.x.npc = 0.4,size=5) +
    theme(plot.title=element_text(size=20,hjust=0.5),
        panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = "", 
        axis.title.x = element_text(size=14),
        axis.title.y = element_text( size=14),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(color = "black")
    )
ggsave(filename = str_glue("{od}/Figure_Group_TIDE_boxplot.pdf"),width=5,height = 5)

#---------------------------IPS analysis------------------------------
ips_df = immune_score(exp=tumor_exprs,method='ips',od=od)$ips
ips_df[1:3,1:3]
plot_data = merge(ips_df,as_tibble(train_model_res$Group,rownames = "sample"),by="sample")
colnames(plot_data)
ggboxplot(plot_data,palette=color_fun1,
    x = "Group", y = "IPS", fill = "Group",
    xlab = "", ylab = "IPS Score",
    title = "", add = ""
) + ylim(7,12)+
    stat_compare_means(label.y.npc = 1, label.x.npc = 0.4,size=5) +
    theme(plot.title=element_text(size=20,hjust=0.5),
        panel.grid = element_blank(), legend.title = element_blank(),
        legend.position = "", 
        axis.title.x = element_text(size=14),
        axis.title.y = element_text( size=14),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.text = element_text(color = "black")
    )
ggsave(filename = str_glue("{od}/Figure_IPS_boxplot.pdf"),width=5,height = 5)




