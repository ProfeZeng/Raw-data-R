using(survival,survminer)
load("data/train_data.RData")

cc <- train_data$data_clinical %>% select(bcr_patient_barcode,OS.Time,OS.Status)
score <- fread("results/3.WGCNA/ssgsea_score.csv") %>% 
    mutate(Sample=str_sub(Sample,1,12))

ss <- merge(cc,score,by.x="bcr_patient_barcode",by.y="Sample")
ss <- mutate(ss,OS.Time=OS.Time/365)
# ss %<>% mutate(Group=ifelse(Score>median(Score),"High","Low"))


res.cut <- surv_cutpoint(ss, #数据集
                         time = "OS.Time", #生存状态
                         event = "OS.Status", #生存时间
                         variables = c("Score") #需要计算的数据列名
                         )
res.cat <- surv_categorize(res.cut)

fit <- survfit(Surv(OS.Time,OS.Status) ~ Score,  # 创建生存对象 
               data = res.cat) # 数据集来源

p <- ggsurvplot(fit, # 创建的拟合对象
           data = ss,  # 指定变量数据来源
        #    conf.int = TRUE, # 显示置信区间
           pval = TRUE, # 添加P值
        #    surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE, # 添加风险表
           xlab = "Years", # 指定x轴标签
           legend = c(0.8,0.75), # 指定图例位置
           legend.title = "Score", # 设置图例标题
           legend.labs = c("High", "Low"), # 指定图例分组标签
           break.x.by = 2
           )  # 设置x轴刻度间距

pdf("results/4.survial/survial.pdf")
print(p,newpage=F)
dev.off()

ggsave(plot=p,"results/4.survial/survial.pdf")
