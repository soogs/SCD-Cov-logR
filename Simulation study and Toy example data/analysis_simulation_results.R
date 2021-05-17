# results compiling and analysis #
# initiated: 17-may-2021 #

# This script compiles the results from the simulation study, and takes a look at some plots

# 1. load data ####
setwd("C:\\Users\\park\\Desktop\\project_binary\\industry\\simulation study\\simulation 11-may-2021-blade\\")
# PATH TO BE CHANGED 
  
result1 <- read.table("./results1_11_may_2021.txt", sep = ",")
result2 <- read.table("./results2_11_may_2021.txt", sep = ",")
result3 <- read.table("./results3_11_may_2021.txt", sep = ",")
result4 <- read.table("./results4_11_may_2021.txt", sep = ",")
result5 <- read.table("./results1_sup_11_may_2021.txt", sep = ",")
result6 <- read.table("./results2_sup_11_may_2021.txt", sep = ",")
result7 <- read.table("./results3_sup_11_may_2021.txt", sep = ",")
result8 <- read.table("./results4_sup_11_may_2021.txt", sep = ",")
result9 <- read.table("./results5_sup_11_may_2021.txt", sep = ",")
result10 <- read.table("./results6_sup_11_may_2021.txt", sep = ",")


# combining the results #
# results1 to results4 go together, while results5 to results10 go together.
# so they have to be separately combined and joined later.

# results1 to results4 #
results <- result1

results[complete.cases(result2),] <- result2[complete.cases(result2),]
results[complete.cases(result3),] <- result3[complete.cases(result3),]
results[complete.cases(result4),] <- result4[complete.cases(result4),]

# results5 to results10 #
results_sup <- result5

results_sup[complete.cases(result6),] <- result6[complete.cases(result6),]
results_sup[complete.cases(result7),] <- result7[complete.cases(result7),]
results_sup[complete.cases(result8),] <- result8[complete.cases(result8),]
results_sup[complete.cases(result9),] <- result9[complete.cases(result9),]
results_sup[complete.cases(result10),] <- result10[complete.cases(result10),]

# checking for NA #
anyNA(results)
anyNA(results_sup)


table(results$signal_level)
table(results_sup$signal_level)


# rbind results and results_sup #
results_final <- rbind(results, results_sup)

dim(results_final)

anyNA(results_final)

table(results_final$dimensions)

table(results_final$py_pattern)

table(results_final$signal_level)

# redefining the results object #
results <- results_final

colnames(results)

tail(results)

# 2. plot making ####

# there are two types of plots: plot1 and plot2.
# plot1 is a general plot where the results from each condition are combined.
# plot2 shows boxplots from each individual condition. 
# In the paper, only plot2 was reported; it contains more information.

ber_df <- results

VAF_factor <- factor(paste("VAF ", ber_df$signal_level, sep = ""), levels = c("VAF 0.8", "VAF 0.5", "VAF 0.2"))

dimension_factor <- factor(paste(ber_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))

relevs <- ber_df$py_pattern
relevs[ber_df$py_pattern == 1] <- "D1, D2"
relevs[ber_df$py_pattern == 2] <- "D1, C"

relevant_factor <- factor(paste("Relevant: ", relevs, sep=""), levels = c("Relevant: D1, D2", "Relevant: D1, C"))

ber_df$dimensions <- dimension_factor2
ber_df$py_pattern <- relevant_factor
ber_df$signal_level <- VAF_factor

# ber_df <- subset(ber_df, sd_x != "sd 50")

which(colnames(ber_df) == "scd_pred_ber")
which(colnames(ber_df) == "diacon_pred_ber")
which(colnames(ber_df) == "diathree_pred_ber")

which(colnames(ber_df) == "diaweighted_pred_ber")
which(colnames(ber_df) == "diaweighted_pred_ber")

which(colnames(ber_df) == "diacor_pred_ber")
which(colnames(ber_df) == "diacor_pred_ber")


scd_ber_index <- which(colnames(ber_df) == "scd_pred_ber")

scd2_ber_index <- which(colnames(ber_df) == "scd2_pred_ber")

diacon_ber_index <- which(colnames(ber_df) == "diacon_pred_ber")

diacor_ber_index <- which(colnames(ber_df) == "diacor_pred_ber")

colnames(ber_df)[c(which(colnames(ber_df) == "scd_pred_ber"),
                   which(colnames(ber_df) == "scd2_pred_ber"),
                   which(colnames(ber_df) == "diacon_pred_ber"), 
                   which(colnames(ber_df) == "diacor_pred_ber"))] <- c("SCD-Cov-logR", "SCD-Cov-logR2" ,"DIABLO-conc", "DIABLO")



ber_df <- ber_df[,c(which(colnames(ber_df) == "dimensions"), which(colnames(ber_df) == "signal_level"), which(colnames(ber_df) == "py_pattern"), scd_ber_index, scd2_ber_index, diacon_ber_index, diacor_ber_index)]

ber_long <- tidyr::gather(data = ber_df, method, ber, c("SCD-Cov-logR", "DIABLO")) # here is where I exclude DIABLO-conc



ber_long_method <- factor(ber_long$method, levels = c("SCD-Cov-logR", "DIABLO"))

ber_long$method <- ber_long_method

library(ggplot2)

# plot 1: general #
plot1_ber <- ggplot(ber_long, aes(x = method, y = ber, fill = method)) +
  theme_bw() +
  labs(fill = "Method") + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  # scale_y_continuous(breaks=seq(0, 1.6, 0.2)) +
  geom_boxplot(width = 0.6, fatten = 1.3, outlier.size = 0.1, lwd = 0.4) +
  stat_summary(fun.y = mean, geom = "point", size = 0.3, shape = 20, color = "red", position = position_dodge(0.6))


# plot 2: conditions #
plot2_ber <- ggplot(ber_long, aes(x = method, y = ber, fill = method)) +
  theme_bw() +
  facet_grid(py_pattern ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) +
  stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6))



# plots for correct classification #
corrects_df <- results

corrects_index <- grep("correct", colnames(corrects_df))

corrects_index <- corrects_index[-4]  # excluding DIABLO-NULL here

# corrects_index <- corrects_index[c(1,3,5)]

colnames(corrects_df)[corrects_index]

corrects_df$dimensions <- dimension_factor2
corrects_df$py_pattern <- relevant_factor
corrects_df$signal_level <- VAF_factor

colnames(corrects_df)[corrects_index] <- c("SCD-Cov-logR", "SCD-Cov-logR2", "DIABLO-conc", "DIABLO-weighted", "DIABLO")

corrects_df <- corrects_df[,c(which(colnames(corrects_df) == "dimensions"), which(colnames(corrects_df) == "signal_level"), which(colnames(corrects_df) == "py_pattern"), corrects_index)]

corrects_long <- tidyr::gather(data = corrects_df, method, corrects, c("SCD-Cov-logR", "DIABLO"))

corrects_long_method <- factor(corrects_long$method, levels = c("SCD-Cov-logR", "DIABLO"))

corrects_long$method <- corrects_long_method

# plot 1: general #
plot1_corrects <- ggplot(corrects_long, aes(x = method, y = corrects, fill = method)) +
  theme_bw() +
  labs(fill = "Method") + 
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  # scale_y_continuous(breaks=seq(0, 1.6, 0.2)) +
  geom_boxplot(width = 0.6, fatten = 1.3, outlier.size = 0.1, lwd = 0.4) +
  stat_summary(fun.y = mean, geom = "point", size = 0.3, shape = 20, color = "red", position = position_dodge(0.6))


# plot 2: conditions #
plot2_corrects <- ggplot(corrects_long, aes(x = method, y = corrects, fill = method)) +
  theme_bw() +
  facet_grid(py_pattern ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title =  element_blank(),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  scale_y_continuous(limits = c(0.65, 1)) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) +
  stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6))


# save(results, file = "C:\\Users\\park\\Desktop\\project_binary\\industry\\simulation study\\simulation 11-may-2021-blade\\results_final_11_may_2021.Rdata")
