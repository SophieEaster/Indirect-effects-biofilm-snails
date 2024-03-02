# Chlorophyll and ash free dry weight of biofilm after colonisation + chronic contamination of 14 days

#packages
library(stats)
library(rstatix)
library(dplyr)
library(tidyverse)
library(tidyr)
library(broom)
library(dunn.test)
library(patchwork)

# chlorophyll -------------------------------------------------------------

# 1) Chlorophyll
# insert data
chlorophyll_biofilm <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\chlorophyll_biofilm.txt',
             header = TRUE, sep = "\t")


View(chlorophyll_biofilm)

# descriptive statistics
#chla
a <- ggplot(chlorophyll_biofilm, aes(x = chla)) +
  geom_histogram(bins = 36)
#chlb
b <- ggplot(chlorophyll_biofilm, aes(x = chlb)) +
  geom_histogram(bins = 36)
#chlc
c <- ggplot(chlorophyll_biofilm, aes(x = chlc)) +
  geom_histogram(bins = 36)

a | b | c

# mean of each treatment
#chl a
mean_chlorophyll_a <- chlorophyll_biofilm %>%
  group_by(week, treatment) %>%
  summarize(mean_chlorophyll_a = mean(chla))

print(mean_chlorophyll_a)

# chl b
mean_chlorophyll_b <- chlorophyll_biofilm %>%
  group_by(week, treatment) %>%
  summarize(mean_chlorophyll_b = mean(chlb))

print(mean_chlorophyll_b)

# chl c
mean_chlorophyll_c <- chlorophyll_biofilm %>%
  group_by(week, treatment) %>%
  summarize(mean_chlorophyll_c = mean(chlc))

print(mean_chlorophyll_c)

# separate table in weeks
week1_chl <-
  chlorophyll_biofilm %>%
  filter(week == "1")
View(week1_chl)

week2_chl <-
  chlorophyll_biofilm %>%
  filter(week == "2")
View(week2_chl)

week3_chl <-
  chlorophyll_biofilm %>%
  filter(week == "3")
View(week3_chl)

# statistical differences within weeks

shapiro_test_chlorophyll <- chlorophyll_biofilm %>%
  group_by(week) %>%
  reframe(
    chla = list(
      W = shapiro.test(chla)$statistic,
      p.value = shapiro.test(chla)$p.value
    ),
    chlb = list(
      W = shapiro.test(chlb)$statistic,
      p.value = shapiro.test(chlb)$p.value
    ),
    chlc = list(
      W = shapiro.test(chlc)$statistic,
      p.value = shapiro.test(chlc)$p.value
    )
  )

View(shapiro_test_chlorophyll)

# all but chl c Week 2 normally distributed
# homoscedasticity test

# week 1 chl a
week1_chl %>%
  summarize(
    statistic = bartlett.test(chla, treatment)$statistic,
    p.value = bartlett.test(chla, treatment)$p.value
  )
# week 1 chl b
week1_chl %>%
  summarize(
    statistic = bartlett.test(chlb, treatment)$statistic,
    p.value = bartlett.test(chlb, treatment)$p.value
  )
# week 1 chl c
week1_chl %>%
  summarize(
    statistic = bartlett.test(chlc, treatment)$statistic,
    p.value = bartlett.test(chlc, treatment)$p.value
  )


# week 2 chla
week2_chl %>%
  summarize(
    statistic = bartlett.test(chla, treatment)$statistic,
    p.value = bartlett.test(chla, treatment)$p.value
  ) # non parametric test
# week 2 chlb
week2_chl %>%
  summarize(
    statistic = bartlett.test(chlb, treatment)$statistic,
    p.value = bartlett.test(chlb, treatment)$p.value
  )


# Week 3 chla
week3_chl %>%
  summarize(
    statistic = bartlett.test(chla, treatment)$statistic,
    p.value = bartlett.test(chla, treatment)$p.value
  )
# Week 3 chlb
week3_chl %>%
  summarize(
    statistic = bartlett.test(chlb, treatment)$statistic,
    p.value = bartlett.test(chlb, treatment)$p.value
  )
# Week 3 chlc
week3_chl %>%
  summarize(
    statistic = bartlett.test(chlc, treatment)$statistic,
    p.value = bartlett.test(chlc, treatment)$p.value
  )




# Anova / kruskal wallis test

# week 1 chla
anova_week1_chla <- aov(week1_chl$chla ~ week1_chl$treatment)
summary(anova_week1_chla)
# no significance

# week 1 chlb
anova_week1_chlb <- aov(week1_chl$chlb ~ week1_chl$treatment)
summary(anova_week1_chlb)
# no significance

# week 1 chlc
anova_week1_chlc <- aov(week1_chl$chlc ~ week1_chl$treatment)
summary(anova_week1_chlc)
# post hoc test
TukeyHSD(aov(anova_week1_chlc))
# no significance


# week 2 chla
week2_chl %>%
  summarize(
    kruskal_statistic = kruskal.test(chla, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(chla, g = treatment)$p.value
  )

# week 2 chlb
anova_week2_chlb <- aov(week2_chl$chlb ~ week2_chl$treatment)
summary(anova_week2_chlb)

# week 2 chl c
week2_chl %>%
  summarize(
    kruskal_statistic = kruskal.test(chlc, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(chlc, g = treatment)$p.value
  )


# week 3 chla
anova_week3_chla <- aov(week3_chl$chla ~ week3_chl$treatment)
summary(anova_week3_chla)

# week 3 chlb
anova_week3_chlb <- aov(week3_chl$chlb ~ week3_chl$treatment)
summary(anova_week3_chlb)

# week 3 chlc
anova_week3_chlc <- aov(week3_chl$chlc ~ week3_chl$treatment)
summary(anova_week3_chlc)


# plots for chlorophyll
# chlorophyll a
ggplot(data = chlorophyll_biofilm, mapping = aes(x = factor(week), y = chla, fill = treatment)) +
  geom_boxplot() +
  scale_fill_okabe_ito() + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = "Chlorophyll a") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

# chlorophyll b
ggplot(data = chlorophyll_biofilm, mapping = aes(x = factor(week), y = chlb, fill = treatment)) +
  geom_boxplot() +
  scale_fill_okabe_ito() + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = "Chlorophyll b") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

# chlorophyll c
ggplot(data = chlorophyll_biofilm, mapping = aes(x = factor(week), y = chlc, fill = treatment)) +
  geom_boxplot() +
  scale_fill_okabe_ito() + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = "Chlorophyll c") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

# AFDW --------------------------------------------------------------------

# insert data
afdw_biofilm_afer_colonisation <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\afdw_after_colonisation.txt',
             header = TRUE, sep = "\t")


View(afdw_biofilm_afer_colonisation)

#descriptive statistics
ggplot(afdw_biofilm_afer_colonisation, aes(x = afdw_after_colonisation)) +
  geom_histogram(bins = 36)

# mean 
mean_afdw_after_colonisation <- afdw_biofilm_afer_colonisation %>%
  group_by(week, treatment) %>%
  summarize(mean_afdw_after_colonisation = mean(afdw_after_colonisation))

print(mean_afdw_after_colonisation)

# separate table in weeks
week1_afdw_after_colonisation <-
  afdw_biofilm_afer_colonisation %>%
  filter(week == "1")
View(week1_afdw_after_colonisation)

week2_afdw_after_colonisation <-
  afdw_biofilm_afer_colonisation %>%
  filter(week == "2")
View(week2_afdw_after_colonisation)

week3_afdw_after_colonisation <-
  afdw_biofilm_afer_colonisation %>%
  filter(week == "3")
View(week3_afdw_after_colonisation)

# statistical differences within weeks
# week 1

week1_afdw_after_colonisation %>%
  summarize(W = shapiro.test(afdw_after_colonisation)$statistic,
            p.value = shapiro.test(afdw_after_colonisation)$p.value)

week1_afdw_after_colonisation %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw_after_colonisation, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw_after_colonisation, g = treatment)$p.value
  )


# week 2
week2_afdw_after_colonisation %>%
  summarize(W = shapiro.test(afdw_after_colonisation)$statistic,
            p.value = shapiro.test(afdw_after_colonisation)$p.value)

week2_afdw_after_colonisation %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw_after_colonisation, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw_after_colonisation, g = treatment)$p.value
  )


# week 3
week3_afdw_after_colonisation %>%
  summarize(W = shapiro.test(afdw_after_colonisation)$statistic,
            p.value = shapiro.test(afdw_after_colonisation)$p.value)

week3_afdw_after_colonisation %>%
  summarize(
    statistic = bartlett.test(afdw_after_colonisation, treatment)$statistic,
    p.value = bartlett.test(afdw_after_colonisation, treatment)$p.value
  )

anova_week3_afdw_after_colonisation <-
  aov(
    week3_afdw_after_colonisation$afdw_after_colonisation ~ week3_afdw_after_colonisation$treatment
  )
summary(anova_week3_afdw_after_colonisation)

TukeyHSD(aov(anova_week3_afdw_after_colonisation))

################### How to analyse between weeks?