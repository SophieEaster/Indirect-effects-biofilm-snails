# feeding rate

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
feeding_rate <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\feeding_rate.txt',
             header = TRUE, sep = "\t")

View(feeding_rate)

feeding_rate <- na.omit(feeding_rate)

View(feeding_rate)

# descriptive statistic
feeding_rate %>%
  get_summary_stats(feeding_rate_ug_day)

ggplot(feeding_rate, aes(x = feeding_rate_ug_day)) +
  geom_histogram(binwidth = 1)

# statistical analyses ----------------------------------------------------

# mean of each treatment
mean_feeding_rate <- feeding_rate %>%
  group_by(week, treatment) %>%
  summarize(mean_feeding_rate = mean(feeding_rate_ug_day))

print(mean_feeding_rate)

# separate table in weeks
week1_feeding_rate <-
  feeding_rate %>%
  filter(week == "1")
View(week1_feeding_rate)

week2_feeding_rate <-
  feeding_rate %>%
  filter(week == "2")
View(week2_feeding_rate)

week3_feeding_rate <-
  feeding_rate %>%
  filter(week == "3")
View(week3_feeding_rate)

# statistical differences within weeks
# week 1 / feeding after 7 days
ggplot(week1_feeding_rate, aes(x = feeding_rate_ug_day)) +
  geom_histogram(binwidth = 1)


week1_feeding_rate %>%
  summarize(W = shapiro.test(feeding_rate_ug_day)$statistic,
            p.value = shapiro.test(feeding_rate_ug_day)$p.value)

week1_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(feeding_rate_ug_day, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(feeding_rate_ug_day, g = treatment)$p.value
  )

week1_feeding_rate %>%
  dunn_test(feeding_rate_ug_day ~ treatment, p.adjust.method = "holm")

# week 2 / feeding after 14 days
ggplot(week2_feeding_rate, aes(x = feeding_rate_ug_day)) +
  geom_histogram(binwidth = 1)

week2_feeding_rate %>%
  summarize(W = shapiro.test(feeding_rate_ug_day)$statistic,
            p.value = shapiro.test(feeding_rate_ug_day)$p.value)

week2_feeding_rate %>%
  summarize(
    statistic = bartlett.test(feeding_rate_ug_day, treatment)$statistic,
    p.value = bartlett.test(feeding_rate_ug_day, treatment)$p.value
  )

week2_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(feeding_rate_ug_day, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(feeding_rate_ug_day, g = treatment)$p.value
  )

week2_feeding_rate %>%
  dunn_test(feeding_rate_ug_day ~ treatment, p.adjust.method = "holm")

# week 3 / feeding after 21 days
ggplot(week3_feeding_rate, aes(x = feeding_rate_ug_day)) +
  geom_histogram(binwidth = 1)

week3_feeding_rate %>%
  summarize(W = shapiro.test(feeding_rate_ug_day)$statistic,
            p.value = shapiro.test(feeding_rate_ug_day)$p.value)

week3_feeding_rate %>%
  summarize(
    statistic = bartlett.test(feeding_rate_ug_day, treatment)$statistic,
    p.value = bartlett.test(feeding_rate_ug_day, treatment)$p.value
  )

week3_feeding_rate %>%
  summarize(
    kruskal_statistic = kruskal.test(feeding_rate_ug_day, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(feeding_rate_ug_day, g = treatment)$p.value
  )

week3_feeding_rate %>%
  dunn_test(feeding_rate_ug_day ~ treatment, p.adjust.method = "holm")

# plots -------------------------------------------------------------------

# boxplot feeding rate after 7 / 14 / 21 days

ggplot(data = feeding_rate, mapping = aes(x = factor(week), y = feeding_rate_ug_day, fill = treatment)) +
  geom_boxplot() +
  scale_fill_okabe_ito() + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = "Feeding Rate") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))
