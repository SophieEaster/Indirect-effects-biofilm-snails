# Fatty acid analyses for snails and NMDS

# load packages
library(vegan)
library(data.table)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(dplyr)
library(stringr)
library(tidyr)
library(shades)
library(ggpubr)
library(base)
library(ggokabeito)
library(gridExtra)
library(multcomp)
library(car)
library(data.table)
library(asbio)
library(plotrix)
library(broom)
library(lmtest)
library(boot)

# FA_Snails_absolute ----------------------------------------------------

# insert data
FA_Snails_absolute <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\FA_Snails_absolute.txt',
             header = TRUE, sep = "\t")


View(FA_Snails_absolute)

# remove NAs, measurement < LOQ
FA_Snails_absolute <- na.omit(FA_Snails_absolute)

# outlier removal, run function one time, after that you can use the shortcut outlierKD (your subset, variable)
outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var), eval(dt))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow = c(2, 2), oma = c(0, 0, 3, 0))
  boxplot(var_name, main = "With outliers")
  hist(var_name, main = "With outliers", xlab = NA, ylab = NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main = "Without outliers")
  hist(var_name, main = "Without outliers", xlab = NA, ylab = NA)
  title("Outlier Check", outer = TRUE)
  na2 <- sum(is.na(var_name))
  cat("Outliers identified:", na2 - na1, "n")
  cat("Propotion (%) of outliers:", round((na2 - na1) /
                                            sum(!is.na(var_name)) * 100, 1), "n")
  cat("Mean of the outliers:", round(mo, 2), "n")
  m2 <- mean(var_name, na.rm = T)
  cat("Mean without removing outliers:", round(m1, 2), "n")
  cat("Mean if we remove outliers:", round(m2, 2), "n")
  
  dt[as.character(substitute(var))] <- invisible(var_name)
  assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
  cat("Outliers successfully removed", "n")
  
}


# single FA comparisons, snail absolute -----------------------------------


# create subgroups to differ between FAs
# C08:0
C08_0 <- subset(FA_Snails_absolute, fa == "C08:0", select = c(replicate, fa, treatment, fa_conc))
View(C08_0)

# Shortcut for outlier
outlierKD(C08_0, fa_conc)
nooutlier_C08_0 <- na.omit(C08_0)
nooutlier_C08_0 
View(nooutlier_C08_0)

# mean
mean_C08_0 <- nooutlier_C08_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C08_0)

# statistics
shapiro.test(nooutlier_C08_0$fa_conc) # not normal distributed
kruskal.test(nooutlier_C08_0$fa_conc ~ nooutlier_C08_0$treatment) # p = 0.3304

# C10:0
# Not considered, too less values.

# C11:0
# Not considered.

# C12:0
C12_0 <- subset(FA_Snails_absolute, fa == "C12:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C12_0, fa_conc)
nooutlier_C12_0 <- na.omit(C12_0)

# mean
mean_C12_0 <- nooutlier_C12_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C12_0)

# statistics
shapiro.test(nooutlier_C12_0$fa_conc) # not normal distributed
kruskal.test(nooutlier_C12_0$fa_conc ~ nooutlier_C12_0$treatment) # p = 0.3586

# C13:0
C13_0 <- subset(FA_Snails_absolute, fa == "C13:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C13_0, fa_conc)
nooutlier_C13_0 <- na.omit(C13_0)

# mean
mean_C13_0 <- nooutlier_C13_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C13_0)

# statistics
shapiro.test(nooutlier_C13_0$fa_conc) # not normal distributed
kruskal.test(nooutlier_C13_0$fa_conc ~ nooutlier_C13_0$treatment) # p = 0.1182

# C14:0
C14_0 <- subset(FA_Snails_absolute, fa == "C14:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C14_0, fa_conc)
nooutlier_C14_0 <- na.omit(C14_0)

# mean
mean_C14_0 <- nooutlier_C14_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C14_0)

# statistics
shapiro.test(nooutlier_C14_0$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C14_0) # variances not equal
kruskal.test(nooutlier_C14_0$fa_conc ~ nooutlier_C14_0$treatment) # p = 0.1479

# C14:1 not considered

# C15:0
C15_0 <- subset(FA_Snails_absolute, fa == "C15:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C15_0, fa_conc)
nooutlier_C15_0 <- na.omit(C15_0)

# mean
mean_C15_0 <- nooutlier_C15_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C15_0)

# statistics
shapiro.test(nooutlier_C15_0$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C15_0) # variances equal
anova_C15_0 <- aov(nooutlier_C15_0$fa_conc~nooutlier_C15_0$treatment)
summary(anova_C15_0)
TukeyHSD(aov(nooutlier_C15_0$fa_conc~nooutlier_C15_0$treatment)) # Mix-Cip: p = 0.02

# C15:1 not considered

# C16:0
C16_0 <- subset(FA_Snails_absolute, fa == "C16:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C16_0, fa_conc)
nooutlier_C16_0 <- na.omit(C16_0)

# mean
mean_C16_0 <- nooutlier_C16_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C16_0)

# statistics
shapiro.test(nooutlier_C16_0$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C16_0) # variances equal
anova_C16_0 <- aov(nooutlier_C16_0$fa_conc~nooutlier_C16_0$treatment)
summary(anova_C16_0)
TukeyHSD(aov(nooutlier_C16_0$fa_conc~nooutlier_C16_0$treatment)) # no significance

# C16:1
C16_1 <- subset(FA_Snails_absolute, fa == "C16:1", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C16_1, fa_conc)
nooutlier_C16_1 <- na.omit(C16_1)

# mean
mean_C16_1 <- nooutlier_C16_1 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C16_1)

# statistics
shapiro.test(nooutlier_C16_1$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C16_1) # variances equal
anova_C16_1 <- aov(nooutlier_C16_1$fa_conc~nooutlier_C16_1$treatment)
summary(anova_C16_1)
TukeyHSD(aov(nooutlier_C16_1$fa_conc~nooutlier_C16_1$treatment)) # NC-Mix: p = 0.0588

# C17:0
C17_0 <- subset(FA_Snails_absolute, fa == "C17:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C17_0, fa_conc)
nooutlier_C17_0 <- na.omit(C17_0)

# mean
mean_C17_0 <- nooutlier_C17_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C17_0)

# statistics
shapiro.test(nooutlier_C17_0$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C17_0) # variances not equal
kruskal.test(nooutlier_C17_0$fa_conc ~ nooutlier_C17_0$treatment) # p = 0.04225
dunn_test(nooutlier_C17_0, fa_conc ~ treatment, p.adjust.method = "holm") #no significance

# C17:1 not considered

# C18:0
C18_0 <- subset(FA_Snails_absolute, fa == "C18:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_0, fa_conc)
nooutlier_C18_0 <- na.omit(C18_0)

# mean
mean_C18_0 <- nooutlier_C18_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_0)

# statistics
shapiro.test(nooutlier_C18_0$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C18_0) # variances not equal
kruskal.test(nooutlier_C18_0$fa_conc ~ nooutlier_C18_0$treatment) 
dunn_test(nooutlier_C18_0, fa_conc ~ treatment, p.adjust.method = "holm") # NC-Cip: p = 0.05859

# C18:1
C18_1 <- subset(FA_Snails_absolute, fa == "C18:1", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_1, fa_conc)
nooutlier_C18_1 <- na.omit(C18_1)

# mean
mean_C18_1 <- nooutlier_C18_1 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_1)

# statistics
shapiro.test(nooutlier_C18_1$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C18_1) # variances equal
anova_C18_1 <- aov(nooutlier_C18_1$fa_conc~nooutlier_C18_1$treatment)
summary(anova_C18_1)
TukeyHSD(aov(nooutlier_C18_1$fa_conc~nooutlier_C18_1$treatment)) # no significance

# C18:2(n-6c)
C18_2c <- subset(FA_Snails_absolute, fa == "C18:2(n-6c)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_2c, fa_conc)
nooutlier_C18_2c <- na.omit(C18_2c)

# mean
mean_C18_2c <- nooutlier_C18_2c %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_2c)

# statistics
shapiro.test(nooutlier_C18_2c$fa_conc) # normal distribution
bartlett.test(fa_conc ~ treatment, nooutlier_C18_2c) # variances equal
anova_C18_2c <- aov(nooutlier_C18_2c$fa_conc~nooutlier_C18_2c$treatment)
summary(anova_C18_2c)
TukeyHSD(aov(nooutlier_C18_2c$fa_conc~nooutlier_C18_2c$treatment)) # NC-Mix: p = 0.0379

# C18:2(n-6t)
C18_2t <- subset(FA_Snails_absolute, fa == "C18:2(n-6t)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_2t, fa_conc)
nooutlier_C18_2t <- na.omit(C18_2t)

# mean
mean_C18_2t <- nooutlier_C18_2t %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_2t)

# statistics
shapiro.test(nooutlier_C18_2t$fa_conc) # not normal distributed
kruskal.test(nooutlier_C18_2t$fa_conc ~ nooutlier_C18_2t$treatment) 
dunn_test(nooutlier_C18_2t, fa_conc ~ treatment, p.adjust.method = "holm") # no significance

# C18:3(n-3)
C18_33 <- subset(FA_Snails_absolute, fa == "C18:3(n-3)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_33, fa_conc)
nooutlier_C18_33 <- na.omit(C18_33)

# mean
mean_C18_33 <- nooutlier_C18_33 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_33)

# statistics
shapiro.test(nooutlier_C18_33$fa_conc) # not normal distributed
kruskal.test(nooutlier_C18_33$fa_conc ~ nooutlier_C18_33$treatment) 
dunn_test(nooutlier_C18_33, fa_conc ~ treatment, p.adjust.method = "holm") # NC-Mix: p = 0.0124*

# C18:3(n-6)
C18_36 <- subset(FA_Snails_absolute, fa == "C18:3(n-6)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_36, fa_conc)
nooutlier_C18_36 <- na.omit(C18_36)

# mean
mean_C18_36 <- nooutlier_C18_36 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_36)

# statistics
shapiro.test(nooutlier_C18_36$fa_conc) # not normal distributed
kruskal.test(nooutlier_C18_36$fa_conc ~ nooutlier_C18_36$treatment) 
dunn_test(nooutlier_C18_36, fa_conc ~ treatment, p.adjust.method = "holm") # no significance

# C18:4(n-3)
C18_4 <- subset(FA_Snails_absolute, fa == "C18:4(n-3)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C18_4, fa_conc)
nooutlier_C18_4 <- na.omit(C18_4)

# mean
mean_C18_4 <- nooutlier_C18_4 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C18_4)

# statistics
shapiro.test(nooutlier_C18_4$fa_conc) # not normal distributed
kruskal.test(nooutlier_C18_4$fa_conc ~ nooutlier_C18_4$treatment) 
dunn_test(nooutlier_C18_4, fa_conc ~ treatment, p.adjust.method = "holm") # no significance

# C20:0
C20_0 <- subset(FA_Snails_absolute, fa == "C20:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C20_0, fa_conc)
nooutlier_C20_0 <- na.omit(C20_0)

# mean
mean_C20_0 <- nooutlier_C20_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C20_0)

# statistics
shapiro.test(nooutlier_C20_0$fa_conc) # not normal distributed
kruskal.test(nooutlier_C20_0$fa_conc ~ nooutlier_C20_0$treatment) 
dunn_test(nooutlier_C20_0, fa_conc ~ treatment, p.adjust.method = "holm") # no significance

# C20:1
C20_1 <- subset(FA_Snails_absolute, fa == "C20:1", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C20_1, fa_conc)
nooutlier_C20_1 <- na.omit(C20_1)

# mean
mean_C20_1 <- nooutlier_C20_1 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C20_1)

# statistics
shapiro.test(nooutlier_C20_1$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C20_1) # variances equal
anova_C20_1 <- aov(nooutlier_C20_1$fa_conc~nooutlier_C20_1$treatment)
summary(anova_C20_1)
TukeyHSD(aov(nooutlier_C20_1$fa_conc~nooutlier_C20_1$treatment)) # NC-Mix: p = 0.0622

# C20:2(n-6)
C20_2 <- subset(FA_Snails_absolute, fa == "C20:2(n-6)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C20_2, fa_conc)
nooutlier_C20_2 <- na.omit(C20_2)

# mean
mean_C20_2 <- nooutlier_C20_2 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C20_2)

# statistics
shapiro.test(nooutlier_C20_2$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C20_2) # variances equal
anova_C20_2 <- aov(nooutlier_C20_2$fa_conc~nooutlier_C20_2$treatment)
summary(anova_C20_2)
TukeyHSD(aov(nooutlier_C20_2$fa_conc~nooutlier_C20_2$treatment)) 
# NC-Cip: p = 0.01831 *
# NC-Mix: p = 0.00090 **
# NC-Pro: p = 0.05316

# C20:3(n-3), not considered

# C20:3(n-6), not considered

# C20:4(n-6)
C20_4 <- subset(FA_Snails_absolute, fa == "C20:4(n-6)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C20_4, fa_conc)
nooutlier_C20_4 <- na.omit(C20_4)

# mean
mean_C20_4 <- nooutlier_C20_4 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C20_4)

# statistics
shapiro.test(nooutlier_C20_4$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C20_4) # variances equal
anova_C20_4 <- aov(nooutlier_C20_4$fa_conc~nooutlier_C20_4$treatment)
summary(anova_C20_4)
TukeyHSD(aov(nooutlier_C20_4$fa_conc~nooutlier_C20_4$treatment)) # no significance

# C20:5(n-3)
C20_5 <- subset(FA_Snails_absolute, fa == "C20:5(n-3)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C20_5, fa_conc)
nooutlier_C20_5 <- na.omit(C20_5)

# mean
mean_C20_5 <- nooutlier_C20_5 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C20_5)

# statistics
shapiro.test(nooutlier_C20_5$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C20_5) # variances equal
anova_C20_5 <- aov(nooutlier_C20_5$fa_conc~nooutlier_C20_5$treatment)
summary(anova_C20_5)
TukeyHSD(aov(nooutlier_C20_5$fa_conc~nooutlier_C20_5$treatment)) # no significance

# C21:0, not considered

# C22:0, not considered

# C22:1
C22_1 <- subset(FA_Snails_absolute, fa == "C22:1", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C22_1, fa_conc)
nooutlier_C22_1 <- na.omit(C22_1)

# mean
mean_C22_1 <- nooutlier_C22_1 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C22_1)

# statistics
shapiro.test(nooutlier_C22_1$fa_conc) # not normal distributed
kruskal.test(nooutlier_C22_1$fa_conc ~ nooutlier_C22_1$treatment) 
dunn_test(nooutlier_C22_1, fa_conc ~ treatment, p.adjust.method = "holm")
# NC-Cip: p = 0.00107 **
# NC-Mix: p = 0.00211 **
# NC-Pro: p = 0.0673

# C22:2
C22_2 <- subset(FA_Snails_absolute, fa == "C22:2", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C22_2, fa_conc)
nooutlier_C22_2 <- na.omit(C22_2)

# mean
mean_C22_2 <- nooutlier_C22_2 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C22_2)

# statistics
shapiro.test(nooutlier_C22_2$fa_conc) # not normal distributed
kruskal.test(nooutlier_C22_2$fa_conc ~ nooutlier_C22_2$treatment) 
dunn_test(nooutlier_C22_2, fa_conc ~ treatment, p.adjust.method = "holm") # no significance

# C22:5(n-3)
C22_5 <- subset(FA_Snails_absolute, fa == "C22:5(n-3)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C22_5, fa_conc)
nooutlier_C22_5 <- na.omit(C22_5)

# mean
mean_C22_5 <- nooutlier_C22_5 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C22_5)

# statistics
shapiro.test(nooutlier_C22_5$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C22_5) # variances equal
anova_C22_5 <- aov(nooutlier_C22_5$fa_conc~nooutlier_C22_5$treatment)
summary(anova_C22_5)
TukeyHSD(aov(nooutlier_C22_5$fa_conc~nooutlier_C22_5$treatment)) # no significance

# C22:6(n-3)
C22_6 <- subset(FA_Snails_absolute, fa == "C22:6(n-3)", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C22_6, fa_conc)
nooutlier_C22_6 <- na.omit(C22_6)

# mean
mean_C22_6 <- nooutlier_C22_6 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C22_6)

# statistics
shapiro.test(nooutlier_C22_6$fa_conc) # normal distributed
bartlett.test(fa_conc ~ treatment, nooutlier_C22_6) # variances equal
anova_C22_6 <- aov(nooutlier_C22_6$fa_conc~nooutlier_C22_6$treatment)
summary(anova_C22_6)
TukeyHSD(aov(nooutlier_C22_6$fa_conc~nooutlier_C22_6$treatment))
# NC-Mix: p = 0.04265*

# C24:0
C24_0 <- subset(FA_Snails_absolute, fa == "C24:0", select = c(replicate, fa, treatment, fa_conc))

# Shortcut for outlier
outlierKD(C24_0, fa_conc)
nooutlier_C24_0 <- na.omit(C24_0)

# mean
mean_C24_0 <- nooutlier_C24_0 %>%
  group_by(treatment) %>%
  summarize(mean_FA_conc = mean(fa_conc, na.rm = TRUE))
print(mean_C24_0)

# statistics
shapiro.test(nooutlier_C24_0$fa_conc) # not normal distributed
kruskal.test(nooutlier_C24_0$fa_conc ~ nooutlier_C24_0$treatment) 
dunn_test(nooutlier_C24_0, fa_conc ~ treatment, p.adjust.method = "holm")
# no significance


# total SAFA, MUFA, PUFA snail absolute ----------------------------------
# total SAFA, MUFA, PUFA content; sa = "snail absolute"
# total SAFA content
SAFA_sa <- bind_rows(nooutlier_C08_0, nooutlier_C12_0, nooutlier_C13_0, nooutlier_C14_0, nooutlier_C15_0, nooutlier_C16_0,
                     nooutlier_C17_0, nooutlier_C18_0, nooutlier_C20_0, nooutlier_C24_0)

total_SAFA_per_snail <- SAFA_sa %>%
  group_by(replicate, treatment) %>%
  summarize(sum_SAFA = sum(fa_conc))
print(total_SAFA_per_snail)

# mean
mean_SAFA_sa <- total_SAFA_per_snail %>%
  group_by(treatment) %>%
  summarize(mean_SAFA = mean(sum_SAFA, na.rm = TRUE))
print(mean_SAFA_sa)

# statistics
shapiro.test(total_SAFA_per_snail $ sum_SAFA) # normal distributed
bartlett.test(sum_SAFA ~ treatment, total_SAFA_per_snail) # variances equal
anova_SAFA_sa <- aov(total_SAFA_per_snail$sum_SAFA ~ total_SAFA_per_snail$treatment)
summary(anova_SAFA_sa)
TukeyHSD(aov(anova_SAFA_sa)) # no significance

# total MUFA content
MUFA_sa <- bind_rows(nooutlier_C16_1, nooutlier_C18_1, nooutlier_C20_1, nooutlier_C22_1)

total_MUFA_per_snail <- MUFA_sa %>%
  group_by(replicate, treatment) %>%
  summarize(sum_MUFA = sum(fa_conc))
print(total_MUFA_per_snail)

# mean
mean_MUFA_sa <- total_MUFA_per_snail %>%
  group_by(treatment) %>%
  summarize(mean_MUFA = mean(sum_MUFA, na.rm = TRUE))
print(mean_MUFA_sa)

# statistics
shapiro.test(total_MUFA_per_snail $ sum_MUFA) # not normal distributed
kruskal.test(total_MUFA_per_snail$sum_MUFA ~ total_MUFA_per_snail$treatment) 
# no significance

# total PUFA content
PUFA_sa <- bind_rows(nooutlier_C18_2c, nooutlier_C18_2t, nooutlier_C18_33, nooutlier_C18_36, nooutlier_C18_4, nooutlier_C20_2,
                     nooutlier_C20_4, nooutlier_C20_5, nooutlier_C22_2, nooutlier_C22_5, nooutlier_C22_6)

total_PUFA_per_snail <- PUFA_sa %>%
  group_by(replicate, treatment) %>%
  summarize(sum_PUFA = sum(fa_conc))
print(total_PUFA_per_snail)

# mean
mean_PUFA_sa <- total_PUFA_per_snail %>%
  group_by(treatment) %>%
  summarize(mean_PUFA = mean(sum_PUFA, na.rm = TRUE))
print(mean_PUFA_sa)

# statistics
shapiro.test(total_PUFA_per_snail $ sum_PUFA) # normal distributed
bartlett.test(sum_PUFA ~ treatment, total_PUFA_per_snail) # variances equal
anova_PUFA_sa <- aov(total_PUFA_per_snail$sum_PUFA ~ total_PUFA_per_snail$treatment)
summary(anova_PUFA_sa)
TukeyHSD(aov(anova_PUFA_sa)) # no significance


# NMDS_Snails_absolute ----------------------------------------------------
#change to wide format
FA_Snails_absolute_wide <-
  pivot_wider(FA_Snails_absolute,
              names_from = fa,
              values_from = fa_conc)

# square root transformation
FA_Snails_absolute_wide <- sqrt(FA_Snails_absolute_wide[, -c(1, 2)])
View(FA_Snails_absolute_wide)

#set NA to 0 
FA_Snails_absolute_wide[is.na(FA_Snails_absolute_wide)] <- 0
View(FA_Snails_absolute_wide)

#meta MDS
NMDS_Snails_absolute <-
  metaMDS(
    FA_Snails_absolute_wide,
    autotransform = FALSE,
    k = 2 ,
    distance = "bray",
    na.rm = TRUE
  )

NMDS_Snails_absolute$stress # stress = 0.1

#add treatment column
colnames(FA_Snails_absolute_wide) <-
  c(
    "C08:0",
    "C10:0",
    "C11:0",
    "C12:0",
    "C13:0",
    "C14:0",
    "C14:1",
    "C15:0",
    "C15:1",
    "C16:0",
    "C16:1",
    "C17:0",
    "C17:1",
    "C18:0",
    "C18:1",
    "C18:2(n-6t)",
    "C18:2(n-6c)",
    "C18:3(n-6)",
    "C18:3(n-3)",
    "C18:4(n-3)",
    "C20:0",
    "C20:1",
    "C20:2(n-6)",
    "C20:3(n-6)",
    "C20:4(n-6)",
    "C20:5(n-3)",
    "C21:0",
    "C22:1",
    "C22:0",
    "C22:5(n-3)",
    "C22:6(n-3)",
    "C24:0",
    "C24:1"
  )

#environmental factors and vectors
# "sa" = snail absolute
FA_variables.SAFAs_sa <-
  envfit(
    NMDS_Snails_absolute,
    subset(
      FA_Snails_absolute_wide,
      select = c(
        "C08:0",
        "C10:0",
        "C11:0",
        "C12:0",
        "C13:0",
        "C14:0",
        "C15:0",
        "C16:0",
        "C17:0",
        "C18:0",
        "C20:0",
        "C21:0",
        "C22:0",
        "C24:0"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.MUFAs_sa <-
  envfit(
    NMDS_Snails_absolute,
    subset(
      FA_Snails_absolute_wide,
      select = c(
        "C14:1",
        "C15:1",
        "C16:1",
        "C17:1",
        "C18:1",
        "C20:1",
        "C22:1",
        "C24:1"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.PUFAs_sa <-
  envfit(
    NMDS_Snails_absolute,
    subset(
      FA_Snails_absolute_wide,
      select = c(
        "C18:2(n-6t)",
        "C18:2(n-6c)",
        "C18:3(n-6)",
        "C18:3(n-3)",
        "C18:4(n-3)",
        "C20:2(n-6)",
        "C20:3(n-6)",
        "C20:4(n-6)",
        "C20:5(n-3)",
        "C22:5(n-3)",
        "C22:6(n-3)"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

#data frame
data.scores_snail_absolute = as.data.frame(scores(NMDS_Snails_absolute)$sites)

Treatments <-
  FA_Snails_absolute $ treatment[FA_Snails_absolute $ fa == "C16:0"]

data.scores_snail_absolute$treatment = Treatments

# vectors for SAFA, MUFA and PUFA
SAFAs_coord_cont_sa = as.data.frame(scores(FA_variables.SAFAs_sa, "vectors")) * ordiArrowMul(FA_variables.SAFAs_sa)
MUFAs_coord_cont_sa = as.data.frame(scores(FA_variables.MUFAs_sa, "vectors")) * ordiArrowMul(FA_variables.MUFAs_sa)
PUFAs_coord_cont_sa = as.data.frame(scores(FA_variables.PUFAs_sa, "vectors")) * ordiArrowMul(FA_variables.PUFAs_sa)


#plot
absol_Snail <- ggplot(data = data.scores_snail_absolute, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(shape = treatment, color = treatment),
    size = 5,
    alpha = 1,
    stroke = 2,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = c("NC" = 16, "Pro" = 17, "Cip" = 18, "Mix" = 15)) +
  scale_color_manual(values = c("NC" = "#009E73", "Pro" = "#E69F00", "Cip" = "#56B4E9", "Mix" = "#F0E442")) +
  stat_conf_ellipse(
    aes(colour = treatment, fill = treatment),
    alpha = 0.3,
    geom = "polygon",
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("NC" = "#009E73", "Pro" = "#E69F00", "Cip" = "#56B4E9", "Mix" = "#F0E442")) +  
  scale_fill_manual(values = c("NC" = "#009E73", "Pro" = "#E69F00", "Cip" = "#56B4E9", "Mix" = "#F0E442")) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1 * 0.1, yend = NMDS2 * 0.1), 
               data = PUFAs_coord_cont_sa, linewidth =1, alpha = 0.5, colour = "black") +
  geom_text_repel(
    data = PUFAs_coord_cont_sa, 
    aes(x = NMDS1 * 0.1, y = NMDS2 * 0.1, label = row.names(PUFAs_coord_cont_sa)), 
    colour = "grey30",
    size = 3,
    segment.color = "black",  # Color of the segment connecting the label to its point
    segment.size = 0.5,  # Size of the segment connecting the label to its point
    direction = "both",  # Allow labels to repel in both directions
    force = 2.5,  # Strength of the repulsion
    box.padding = unit(0.15, "lines")  # Padding around labels
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 10, colour = "grey30"),
    plot.margin = margin(
      t = 0.2,
      b = 0.2,
      l = 0.2,
      r = 0.2,
      unit = "cm"
    ),
    panel.grid.major = element_line(colour = "grey95", size = 0.4),
    panel.grid.minor = element_line(colour = "grey95", size = 0.2)
  ) #+
  #annotate(
  #  "text",
  #  x = -0.23,
  #  y = -0.19,
  #  label = paste("stress:", NMDS_Snails_absolute$stress %>% round(2)),
  #  size = 3.5
  #)


absol_Snail
# position of FA names manually adjusted in vector graphic

# save the plot
ggsave(absol_Snail, 
       filename = "NMDS_snails_absolute.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')



# FA_Snails_relative ----------------------------------------------------
# calculate the distribution of fatty acids in snails in percent
FA_Snails_relative <- FA_Snails_absolute %>%
  group_by(replicate) %>%
  mutate(fa_conc_rel = (fa_conc / sum(fa_conc)) * 100)
print(FA_Snails_relative)

# now single FA comparisons, snails relative

# noch machen

# NMDS_Snails_relative ----------------------------------------------------
#change to wide format
FA_Snails_relative_wide <-
  pivot_wider(FA_Snails_relative,
              names_from = fa,
              values_from = fa_rel)

# square root transformation
FA_Snails_relative_wide <- sqrt(FA_Snails_relative_wide[, -c(1, 2)])
View(FA_Snails_relative_wide)

#set NA to 0 
FA_Snails_relative_wide[is.na(FA_Snails_relative_wide)] <- 0
View(FA_Snails_relative_wide)
#meta MDS
NMDS_Snails_relative <-
  metaMDS(
    FA_Snails_relative_wide,
    autotransform = FALSE,
    k = 2 ,
    distance = "bray",
    na.rm = TRUE
  )

NMDS_Snails_relative$stress

#add treatment column
colnames(FA_Snails_relative_wide) <-
  c(
    "C08:0",
    "C12:0",
    "C13:0",
    "C14:0",
    "C15:0",
    "C16:0",
    "C16:1",
    "C17:0",
    "C18:0",
    "C18:1",
    "C18:2(n-6t)",
    "C18:2(n-6c)",
    "C18:3(n-6)",
    "C18:3(n-3)",
    "C18:4(n-3)",
    "C20:0",
    "C20:1",
    "C20:2(n-6)",
    "C20:4(n-6)",
    "C20:5(n-3)",
    "C22:1",
    "C22:2",
    "C22:5(n-3)",
    "C22:6(n-3)",
    "C24:0"
  )

#environmental factors and vectors
# "sr" = snail relative
FA_variables.SAFAs_sr <-
  envfit(
    NMDS_Snails_relative,
    subset(
      FA_Snails_relative_wide,
      select = c(
        "C08:0",
        "C12:0",
        "C13:0",
        "C14:0",
        "C15:0",
        "C16:0",
        "C17:0",
        "C18:0",
        "C20:0",
        "C24:0"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.MUFAs_sr <-
  envfit(
    NMDS_Snails_relative,
    subset(
      FA_Snails_relative_wide,
      select = c(
        "C16:1",
        "C18:1",
        "C20:1",
        "C22:1"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.PUFAs_sr <-
  envfit(
    NMDS_Snails_relative,
    subset(
      FA_Snails_relative_wide,
      select = c(
        "C18:2(n-6t)",
        "C18:2(n-6c)",
        "C18:3(n-6)",
        "C18:3(n-3)",
        "C18:4(n-3)",
        "C20:2(n-6)",
        "C20:4(n-6)",
        "C20:5(n-3)",
        "C22:5(n-3)",
        "C22:6(n-3)"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

#data frame
data.scores_snail_relative = as.data.frame(scores(NMDS_Snails_relative)$sites)

Treatments <-
  FA_Snails_relative $ treatment[FA_Snails_relative $ fa == "C16:0"]

data.scores_snail_relative$treatment = Treatments

# vectors for SAFA, MUFA and PUFA
SAFAs_coord_cont_sr = as.data.frame(scores(FA_variables.SAFAs_sr, "vectors")) * ordiArrowMul(FA_variables.SAFAs_sr)
MUFAs_coord_cont_sr = as.data.frame(scores(FA_variables.MUFAs_sr, "vectors")) * ordiArrowMul(FA_variables.MUFAs_sr)
PUFAs_coord_cont_sr = as.data.frame(scores(FA_variables.PUFAs_sr, "vectors")) * ordiArrowMul(FA_variables.PUFAs_sr)


#plot
rel_Snail <- ggplot(data = data.scores_snail_relative, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    aes(shape = treatment, color = treatment),
    size = 5,
    alpha = 1,
    stroke = 2,
    show.legend = TRUE
  ) +
  scale_shape_manual(values = c("NC" = 16, "Pro" = 17, "Cip" = 18, "Mix" = 15)) +
  scale_color_manual(values = c("NC" = "#009E73", "Pro" = "#E69F00", "Cip" = "#56B4E9", "Mix" = "#F0E442")) +
  stat_conf_ellipse(
    aes(colour = treatment, fill = treatment),
    alpha = 0.3,
    geom = "polygon",
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("NC" = "#009E73", "Pro" = "#E69F00", "Cip" = "#56B4E9", "Mix" = "#F0E442")) +  
  scale_fill_manual(values = c("NC" = "#009E73", "Pro" = "#E69F00", "Cip" = "#56B4E9", "Mix" = "#F0E442")) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1 * 0.1, yend = NMDS2 * 0.1), 
               data = PUFAs_coord_cont_sr, linewidth =1, alpha = 0.5, colour = "black") +
  geom_text_repel(
    data = PUFAs_coord_cont_sr, 
    aes(x = NMDS1 * 0.1, y = NMDS2 * 0.1, label = row.names(PUFAs_coord_cont_sr)), 
    colour = "grey30",
    size = 3,
    segment.color = "black",  # Color of the segment connecting the label to its point
    segment.size = 0.5,  # Size of the segment connecting the label to its point
    direction = "both",  # Allow labels to repel in both directions
    force = 2.5,  # Strength of the repulsion
    box.padding = unit(0.15, "lines")  # Padding around labels
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 10, colour = "grey30"),
    plot.margin = margin(
      t = 0.2,
      b = 0.2,
      l = 0.2,
      r = 0.2,
      unit = "cm"
    ),
    panel.grid.major = element_line(colour = "grey95", size = 0.4),
    panel.grid.minor = element_line(colour = "grey95", size = 0.2)
  ) #+
  #annotate(
  # "text",
  #x = -0.05,
  #y =  0.17,
  #label = paste("stress:", NMDS_Snails_relative$stress %>% round(2)),
  #size = 3.5
  #)


rel_Snail
# position of FA names manually adjusted in vector graphic

# save the plot
ggsave(rel_Snail, 
       filename = "NMDS_snails_relative.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')





