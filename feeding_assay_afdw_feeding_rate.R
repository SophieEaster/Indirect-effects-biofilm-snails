# feeding assay: ash free dry weight (afdw) and feeding rate

# packages
library(stats)
library(rstatix)
library(dplyr)
library(tidyverse)
library(tidyr)
library(broom)
library(dunn.test)
library(ggplot2)
library(ggokabeito)
library(moments)
library(lme4)
library(DAAG)
library(performance)
library(lsmeans)

# insert data
afdw_feeding_rate <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\feeding_assay_afdw_feeding_rate.txt',
             header = TRUE, sep = "\t")

View(afdw_feeding_rate)

# NA's due to missing snail 
# Zero values of feeding rate due to no measurable feeding in respective week
# unit feeding rate: µg / day
# unit ash free dry weight (afdw): mg / cm2
afdw_feeding_rate <- na.omit(afdw_feeding_rate)
View(afdw_feeding_rate)

# stats feeding rate ----------------------------------------------------

# mean of each treatment
mean_feeding_rate <- afdw_feeding_rate %>%
  group_by(week, treatment) %>%
  summarize(mean_feeding_rate = mean(feeding_rate))

print(mean_feeding_rate)

# separate table in weeks
week1_feeding_rate <-
  afdw_feeding_rate %>%
  filter(week == "1")
View(week1_feeding_rate)

week2_feeding_rate <-
  afdw_feeding_rate %>%
  filter(week == "2")
View(week2_feeding_rate)

week3_feeding_rate <-
  afdw_feeding_rate %>%
  filter(week == "3")
View(week3_feeding_rate)

# statistical differences within weeks
# week 1 / feeding after 7 days
ggplot(week1_feeding_rate, aes(x = feeding_rate)) +
  geom_histogram(binwidth = 1)

week1_feeding_rate %>%
  summarize(W = shapiro.test(feeding_rate)$statistic,
            p.value = shapiro.test(feeding_rate)$p.value)

week1_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(feeding_rate, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(feeding_rate, g = treatment)$p.value
  )

week1_feeding_rate %>%
  dunn_test(feeding_rate ~ treatment, p.adjust.method = "holm")

# week 2 / feeding after 14 days
ggplot(week2_feeding_rate, aes(x = feeding_rate)) +
  geom_histogram(binwidth = 1)

week2_feeding_rate %>%
  summarize(W = shapiro.test(feeding_rate)$statistic,
            p.value = shapiro.test(feeding_rate)$p.value)

week2_feeding_rate %>%
  summarize(
    statistic = bartlett.test(feeding_rate, treatment)$statistic,
    p.value = bartlett.test(feeding_rate, treatment)$p.value
  )

week2_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(feeding_rate, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(feeding_rate, g = treatment)$p.value
  )

week2_feeding_rate %>%
  dunn_test(feeding_rate ~ treatment, p.adjust.method = "holm")

# week 3 / feeding after 21 days
ggplot(week3_feeding_rate, aes(x = feeding_rate)) +
  geom_histogram(binwidth = 1)

week3_feeding_rate %>%
  summarize(W = shapiro.test(feeding_rate)$statistic,
            p.value = shapiro.test(feeding_rate)$p.value)

week3_feeding_rate %>%
  summarize(
    statistic = bartlett.test(feeding_rate, treatment)$statistic,
    p.value = bartlett.test(feeding_rate, treatment)$p.value
  )

week3_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(feeding_rate, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(feeding_rate, g = treatment)$p.value
  )

week3_feeding_rate %>%
  dunn_test(feeding_rate ~ treatment, p.adjust.method = "holm")

# stats afdw after feeding ------------------------------------------------
# mean of each treatment
mean_afdw_after_feeding <- afdw_feeding_rate %>%
  group_by(week, treatment) %>%
  summarize(mean_afdw_after_feeding = mean(afdw))

print(mean_afdw_after_feeding)

# statistical differences within weeks
# week 1 
week1_feeding_rate %>%
  summarize(W = shapiro.test(afdw)$statistic,
            p.value = shapiro.test(afdw)$p.value)

week1_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw, g = treatment)$p.value
  )

week1_feeding_rate %>%
  dunn_test(afdw ~ treatment, p.adjust.method = "holm")

# week 2
week2_feeding_rate %>%
  summarize(W = shapiro.test(afdw)$statistic,
            p.value = shapiro.test(afdw)$p.value)

week2_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw, g = treatment)$p.value
  )

week2_feeding_rate %>%
  dunn_test(afdw ~ treatment, p.adjust.method = "holm")

# week 3 
week3_feeding_rate %>%
  summarize(W = shapiro.test(afdw)$statistic,
            p.value = shapiro.test(afdw)$p.value)

week3_feeding_rate %>%
  summarize(
    statistic = bartlett.test(afdw, treatment)$statistic,
    p.value = bartlett.test(afdw, treatment)$p.value
  )

anova_week3_afdw_after_feeding <- aov(week3_feeding_rate $ afdw ~ week3_feeding_rate $ treatment)
summary(anova_week3_afdw_after_feeding)

TukeyHSD(aov(week3_feeding_rate $ afdw ~ week3_feeding_rate $ treatment))

# plots feeding rate  -------------------------------------------------------------------
# histogram feeding rate for all weeks
ggplot(afdw_feeding_rate, aes(x = feeding_rate)) +
  geom_histogram(binwidth = 1)

# boxplot feeding rate after 7 / 14 / 21 days
feeding_boxplot <- ggplot(data = afdw_feeding_rate, 
       mapping = aes(x = factor(week), 
                     y = feeding_rate, 
                     fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = "Feeding Rate (µg/day)", fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

feeding_boxplot

# save the plot
ggsave(feeding_boxplot, 
       filename = "feeding_rate_boxplot.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')


# plots afdw --------------------------------------------------------------
# histogram afdw for all weeks
ggplot(afdw_feeding_rate, aes(x = afdw)) +
  geom_histogram(bins = 330)

# boxplot feeding rate after 7 / 14 / 21 days
afdw_after_feeding_boxplot <- ggplot(data = afdw_feeding_rate, 
       mapping = aes(x = factor(week), 
                     y = afdw, 
                     fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = expression("AFDW (mg/cm"^"2"*")"), fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

afdw_after_feeding_boxplot

# save the plot
ggsave(afdw_after_feeding_boxplot, 
       filename = "afdw_after_feeding_boxplot.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')



# LMM feeding rate --------------------------------------------------------
# Linear mixed model for feeding rate 
# response variable: feeding rate, explanatory variables: treatment & timepoint (week)

# timepoint week as factor and prefix "week"
afdw_feeding_rate $ week <- paste("week", afdw_feeding_rate $ week, sep = "")

afdw_feeding_rate $ week <- factor(afdw_feeding_rate $ week)

ggplot(afdw_feeding_rate, aes(x = feeding_rate)) +
  geom_histogram(binwidth = 1)
# due to right skewed distribution, log 10 transformation of "feeding_rate" and 
# first remove all rows in which feeding rate = 0 
afdw_feeding_rate <- afdw_feeding_rate %>%
  filter(feeding_rate != 0) %>%
  mutate(log_feeding_rate = log10(feeding_rate))

ggplot(afdw_feeding_rate, aes(x = log_feeding_rate)) +
  geom_histogram(bins = 10)
# normal distributed
# check skewness
skewness(afdw_feeding_rate $ log_feeding_rate)
# 0.2735669
# QQ Plot
qqnorm(afdw_feeding_rate $ log_feeding_rate, datax = TRUE)
qqline(afdw_feeding_rate $ log_feeding_rate, datax = TRUE)
# not perfect fit, but normal distribution assumed

# model building 
# response variable: log_feeding_rate
# fixed effects: week + treatment + week:treatment
# random effect with constant intercept: (1|replicate)
# grouping factor: replicate as factor
afdw_feeding_rate $ replicate <- factor(afdw_feeding_rate $ replicate)
feeding_rate_lmer <- lmer(log_feeding_rate ~ week + treatment + week:treatment + (1|replicate), data = afdw_feeding_rate)

feeding_rate_lmer %>% r2
# Conditional R2: explained portion of the variance of the fixed and random effect together
# Marginal R2: declared share of the variance of the fixed components alone

plot(feeding_rate_lmer) # point cloud looks evenly distributed --> normal distribution

qreference(residuals(feeding_rate_lmer)) # normality

drop1(feeding_rate_lmer, test = 'Chisq')
# interaction of timepoint (week) and treatment highly significant
# now check at which timepoints is a difference between treatments
anova(feeding_rate_lmer)

# conduct Tukey HSD-test for pairwise comparisons
feeding_rate_tukey <- lsmeans(feeding_rate_lmer, ~ treatment | week)
feeding_rate_tukey

feeding_rate_tukey1 <- contrast(feeding_rate_tukey, 'pairwise')
feeding_rate_tukey1

feeding_rate_summary <- summary(feeding_rate_tukey1, infer = TRUE, adjust = 'mvt')
feeding_rate_summary
