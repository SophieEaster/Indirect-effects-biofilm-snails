# Ash free dry weight of biofilm after 7 days of feeding
#packages
library(stats)
library(rstatix)
library(dplyr)
library(tidyverse)
library(tidyr)
library(broom)
library(dunn.test)
#insert data
afdw_after_feeding <-
  read.table(file = 'C:\\Users\\sophi\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\afdw_after_feeding.txt',
             header = TRUE, sep = "\t")

View(afdw_after_feeding)

afdw_after_feeding <- na.omit(afdw_after_feeding)

View(afdw_after_feeding)

# mean of each treatment
mean_afdw_after_feeding <- afdw_after_feeding %>%
  group_by(week, treatment) %>%
  summarize(mean_afdw_after_feeding = mean(afdw))

print(mean_afdw_after_feeding)

# separate table in weeks
week1_afdw_after_feeding <-
  afdw_after_feeding %>%
  filter(week == "1")
View(week1_afdw_after_feeding)

week2_afdw_after_feeding <-
  afdw_after_feeding %>%
  filter(week == "2")
View(week2_afdw_after_feeding)

week3_afdw_after_feeding <-
  afdw_after_feeding %>%
  filter(week == "3")
View(week3_afdw_after_feeding)

# statistical differences within weeks
# week 1

week1_afdw_after_feeding %>%
  summarize(W = shapiro.test(afdw)$statistic,
            p.value = shapiro.test(afdw)$p.value)

week1_afdw_after_feeding %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw, g = treatment)$p.value
  )

week1_afdw_after_feeding %>%
  dunn_test(afdw ~ treatment, p.adjust.method = "holm")

# week2
week2_afdw_after_feeding %>%
  summarize(W = shapiro.test(afdw)$statistic,
            p.value = shapiro.test(afdw)$p.value)

week2_afdw_after_feeding %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw, g = treatment)$p.value
  )

week2_afdw_after_feeding %>%
  dunn_test(afdw ~ treatment, p.adjust.method = "holm")

# week 3
week3_afdw_after_feeding %>%
  summarize(W = shapiro.test(afdw)$statistic,
            p.value = shapiro.test(afdw)$p.value)

week3_afdw_after_feeding %>%
  summarize(
    statistic = bartlett.test(afdw, treatment)$statistic,
    p.value = bartlett.test(afdw, treatment)$p.value
  )

anova_week3_afdw_after_feeding <- aov(week3_afdw_after_feeding$afdw ~ week3_afdw_after_feeding$treatment)
summary(anova_week3_afdw_after_feeding)

TukeyHSD(aov(week3_afdw_after_feeding$afdw ~ week3_afdw_after_feeding$treatment))

################################################# 
# How to analyse between weeks?
##