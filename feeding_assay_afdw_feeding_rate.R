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
# histogram feeding rate for all weeks
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
