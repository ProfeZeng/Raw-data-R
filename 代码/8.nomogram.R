rm(list=ls())
# 导入函数
source("src/0.config.R")
conflict_prefer("between", "dplyr",quiet = TRUE)
library(survival)
# 设置结果路径
od <- "results/8.nomogram"
suppressWarnings(dir.create(od,recursive=TRUE))

# 导入数据
load("data/train_data.RData")
load("results/5.model/signature/train_model_res.RData")
curves_data <- train_model_res$Score %>% as_tibble(rownames = "Sample")
cc <- train_data$data_clinical %>% select(Sample,OS.Status,OS.Time, Sex,Age,Sex,Grade)
curves_data <- merge(cc,curves_data,by="Sample")

head(curves_data)
# -----------------------------------nomogram---------------------------------
conflict_prefer("filter", "dplyr")
# GEO数据集
nomogram_df <- curves_data[,-1]
head(nomogram_df)

model <- survival::coxph(Surv(OS.Time, OS.Status) ~ Sex+Age+Grade+Score, data=nomogram_df)
library(regplot)
pdf(paste0(od, "/nomogram.pdf"), w = 5, h = 5)
p=regplot(model, plots = c("violin", "boxes"), observation = TRUE, title = "Nomogram", clickable = TRUE, points = TRUE, droplines = F)
dev.off()

# ---------- DCA ---------------
install.packages("ggDCA")
DCA <- curves_data
DCA$Sample <- NULL
DCA$Sex %<>% factor()
DCA$Age %<>% factor()
DCA$Grade %<>% factor()

# remotes::install_github('yikeshu0611/ggDCA')
library(ggDCA)
library(rms)
ddist <- rms::datadist(DCA)  ###打包数据
options(datadist = "ddist")

f1<-rms::cph(Surv(OS.Time, OS.Status)~Sex+Age+Grade+Score,DCA)
dca_cph <- dca(f1, model.names = c("model"),times = "median")
ggplot(dca_cph, lwd = 1)
ggsave(paste0(od, "/DCA.pdf"), width = 8, height = 6)


# -----------------------------------calibration---------------------------------
calibration_df <- curves_data
calibration_df$Sample <- NULL
p_load(riskRegression, penalized, prodlim)
library(survival)
library(rms)
library(prodlim)
set.seed(18)

pdf(paste0(od, "/Calibration.pdf"), w = 6, h = 6)

cph <- rms::cph(Surv(OS.Time, OS.Status)~Sex+Age+Grade+Score, data = calibration_df, x = TRUE,time.inc=365,y = TRUE, surv = TRUE)
cal <- calibrate(cph, cmethod = "KM", method = "boot", u = 365*1, m = 50, B = 200)

plot(cal, lwd = 2, lty = 1, errbar.col = color_fun1[1],xlab = "Nomogram-Predicted Probability", ylab = "Actual OS (proportion)",col = color_fun1[1], subtitles = FALSE, xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot")
lines(cal[, c("mean.predicted", "KM")], type = "l", lwd = 2, col = color_fun1[1], pch = 16,add=T)

cph <- rms::cph(Surv(OS.Time, OS.Status)~ Sex+Age+Grade+Score, data = calibration_df, x = TRUE,time.inc=365*3,y = TRUE, surv = TRUE)
cal <- rms::calibrate(cph, cmethod = "KM", method = "boot", u = 365*3, m = 50, B = 300)

plot(cal, lwd = 2, lty = 1, errbar.col = color_fun1[2],xlab = "Nomogram-Predicted Probability", ylab = "Actual OS (proportion)",col = color_fun1[2], subtitles = FALSE, xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot",add=T)
lines(cal[, c("mean.predicted", "KM")], type = "l", lwd = 2, col = color_fun1[2], pch = 16)

cph <- rms::cph(Surv(OS.Time, OS.Status)~ Sex+Age+Grade+Score, data = calibration_df, x = TRUE,time.inc=365*5,y = TRUE, surv = TRUE)
cal <- rms::calibrate(cph, cmethod = "KM", method = "boot", u = 365*3, m = 50, B = 300)

plot(cal, lwd = 2, lty = 1, errbar.col = color_fun1[3],xlab = "Nomogram-Predicted Probability", ylab = "Actual OS (proportion)",col = color_fun1[3], subtitles = FALSE, xlim = c(0,1), ylim = c(0, 1), main = "Calibrate plot",add=T)
lines(cal[, c("mean.predicted", "KM")], type = "l", lwd = 2, col = color_fun1[3], pch = 16)
abline(0, 1, lty = 3, lwd = 2, col = "grey")
legend(0.8, 0.2, legend = c("1 year", "3 year","5 year"), col = color_fun1[1:3], lty = 1, cex = 0.8, box.lty = 0)

dev.off()

# -----------------------------------ROC---------------------------------
dat <- curves_data[,-1]
model <- survival::coxph(Surv(OS.Time, OS.Status) ~  Sex+Age+Grade+Score, data=dat)
dat$risk <- predict(model,dat,type="risk")
library(timeROC)
roc_model <- timeROC(
    T = dat$OS.Time,
    delta = dat$OS.Status,
    marker = dat$risk,
    cause = 1,
    weighting = "aalen", #marginal\aalen\cox
    times = c(365,365*3,365*5),
        ROC = TRUE,
        iid = FALSE
)
roc_plot <- data.frame(
year1x = roc_model$FP[, 1], year1y = roc_model$TP[, 1],
year2x = roc_model$FP[, 2], year2y = roc_model$TP[, 2],
year3x = roc_model$FP[, 3], year3y = roc_model$TP[, 3]
)
AUC_anno_1 <- sprintf("AUC at %s year = %s", 1, sprintf("%.3f", roc_model$AUC[[1]]))
AUC_anno_2 <- sprintf("AUC at %s year = %s", 3, sprintf("%.3f", roc_model$AUC[[2]]))
AUC_anno_3 <- sprintf("AUC at %s year = %s", 5, sprintf("%.3f", roc_model$AUC[[3]]))

p2 <- ggplot(data = roc_plot) +
geom_line(aes(x = year1x, y = year1y), size = 1.2, color = "#E31A1C") +
geom_line(aes(x = year2x, y = year2y), size = 1.2, color = "#377EB8") +
geom_line(aes(x = year3x, y = year3y), size = 1.2, color = "#007947") +
geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
egg::theme_article() +
theme(plot.background = element_rect(fill = "white")) +
annotate(geom = "line", x = c(0.5, 0.54), y = .17, colour = "#E31A1C", size = 1.2) +
annotate("text", x = 0.55, y = .17, size = 5, label = AUC_anno_1, color = "black", hjust = "left") +
annotate(geom = "line", x = c(0.5, 0.54), y = .11, colour = "#377EB8", size = 1.2) +
annotate("text", x = 0.55, y = .11, size = 5, label = AUC_anno_2, color = "black", hjust = "left") +
annotate(geom = "line", x = c(0.5, 0.54), y = .05, colour = "#007947", size = 1.2) +
annotate("text", x = 0.55, y = .05, size = 5, label = AUC_anno_3, color = "black", hjust = "left") +
labs(x = "1-Specificity", y = "Sensitivity") +
theme(
    axis.text.x = element_text(face = "plain", size = 12, color = "black"),
    axis.text.y = element_text(face = "plain", size = 12, color = "black"),
    axis.title.x = element_text(face = "plain", size = 14, color = "black"),
    axis.title.y = element_text(face = "plain", size = 14, color = "black")
)

ggsave(paste0(od, "/ROC.pdf"), width  = 6, height  = 6)

