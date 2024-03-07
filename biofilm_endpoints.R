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

# insert data
biofilm_endpoints <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\biofilm_endpoints.txt',
             header = TRUE, sep = "\t")

View(biofilm_endpoints)

# Chlorophyll -------------------------------------------------------------

# mean of each treatment
#chl a
mean_chlorophyll_a <- biofilm_endpoints %>%
  group_by(week, treatment) %>%
  summarize(mean_chlorophyll_a = mean(chla))

print(mean_chlorophyll_a)

# chl b
mean_chlorophyll_b <- biofilm_endpoints %>%
  group_by(week, treatment) %>%
  summarize(mean_chlorophyll_b = mean(chlb))

print(mean_chlorophyll_b)

# chl c
mean_chlorophyll_c <- biofilm_endpoints %>%
  group_by(week, treatment) %>%
  summarize(mean_chlorophyll_c = mean(chlc))

print(mean_chlorophyll_c)

# separate table in weeks
biofilm_endpoints_week1 <-
  biofilm_endpoints %>%
  filter(week == "1") 
View(biofilm_endpoints_week1)

biofilm_endpoints_week2 <-
  biofilm_endpoints %>%
  filter(week == "2")
View(biofilm_endpoints_week2)

biofilm_endpoints_week3 <-
  biofilm_endpoints %>%
  filter(week == "3")
View(biofilm_endpoints_week3)

# statistical differences within weeks

shapiro_test_chlorophyll <- biofilm_endpoints %>%
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
biofilm_endpoints_week1 %>%
  summarize(
    statistic = bartlett.test(chla, treatment)$statistic,
    p.value = bartlett.test(chla, treatment)$p.value
  )
# week 1 chl b
biofilm_endpoints_week1 %>%
  summarize(
    statistic = bartlett.test(chlb, treatment)$statistic,
    p.value = bartlett.test(chlb, treatment)$p.value
  )
# week 1 chl c
biofilm_endpoints_week3 %>%
  summarize(
    statistic = bartlett.test(chlc, treatment)$statistic,
    p.value = bartlett.test(chlc, treatment)$p.value
  )


# week 2 chla
biofilm_endpoints_week2 %>%
  summarize(
    statistic = bartlett.test(chla, treatment)$statistic,
    p.value = bartlett.test(chla, treatment)$p.value
  ) # non parametric test
# week 2 chlb
biofilm_endpoints_week2 %>%
  summarize(
    statistic = bartlett.test(chlb, treatment)$statistic,
    p.value = bartlett.test(chlb, treatment)$p.value
  )


# Week 3 chla
biofilm_endpoints_week3 %>%
  summarize(
    statistic = bartlett.test(chla, treatment)$statistic,
    p.value = bartlett.test(chla, treatment)$p.value
  )
# Week 3 chlb
biofilm_endpoints_week3 %>%
  summarize(
    statistic = bartlett.test(chlb, treatment)$statistic,
    p.value = bartlett.test(chlb, treatment)$p.value
  )
# Week 3 chlc
biofilm_endpoints_week3 %>%
  summarize(
    statistic = bartlett.test(chlc, treatment)$statistic,
    p.value = bartlett.test(chlc, treatment)$p.value
  )


# Anova / kruskal wallis test

# week 1 chla
anova_week1_chla <- aov(biofilm_endpoints_week1 $ chla ~ biofilm_endpoints_week1 $ treatment)
summary(anova_week1_chla)
# no significance

# week 1 chlb
anova_week1_chlb <- aov(biofilm_endpoints_week1 $ chlb ~ biofilm_endpoints_week1 $ treatment)
summary(anova_week1_chlb)
# no significance

# week 1 chlc
anova_week1_chlc <- aov(biofilm_endpoints_week1 $ chlc ~ biofilm_endpoints_week1 $ treatment)
summary(anova_week1_chlc)
# post hoc test
TukeyHSD(aov(anova_week1_chlc))
# no significance


# week 2 chla
biofilm_endpoints_week2 %>%
  summarize(
    kruskal_statistic = kruskal.test(chla, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(chla, g = treatment)$p.value
  )

# week 2 chlb
anova_week2_chlb <- aov(biofilm_endpoints_week2 $ chlb ~ biofilm_endpoints_week2 $ treatment)
summary(anova_week2_chlb)

# week 2 chl c
biofilm_endpoints_week2 %>%
  summarize(
    kruskal_statistic = kruskal.test(chlc, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(chlc, g = treatment)$p.value
  )


# week 3 chla
anova_week3_chla <- aov(biofilm_endpoints_week3 $ chla ~ biofilm_endpoints_week3 $ treatment)
summary(anova_week3_chla)

# week 3 chlb
anova_week3_chlb <- aov(biofilm_endpoints_week3 $ chlb ~ biofilm_endpoints_week3 $ treatment)
summary(anova_week3_chlb)

# week 3 chlc
anova_week3_chlc <- aov(biofilm_endpoints_week3 $ chlc ~ biofilm_endpoints_week3 $ treatment)
summary(anova_week3_chlc)


# Plots chlorophyll -------------------------------------------------------
# descriptive statistics
#chla
a <- ggplot(biofilm_endpoints, aes(x = chla)) +
  geom_histogram(bins = 36)
#chlb
b <- ggplot(biofilm_endpoints, aes(x = chlb)) +
  geom_histogram(bins = 36)
#chlc
c <- ggplot(biofilm_endpoints, aes(x = chlc)) +
  geom_histogram(bins = 36)

a | b | c

# chlorophyll a
chla <- ggplot(data = biofilm_endpoints, 
       mapping = aes(x = factor(week), 
                     y = chla, 
                     fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = expression("Chlorophyll a (µg/cm"^"2"*")"), fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

chla 

# save the plot
ggsave(chla, 
       filename = "chla.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')


# chlorophyll b
chlb <- ggplot(data = biofilm_endpoints, 
               mapping = aes(x = factor(week), 
                             y = chlb, 
                             fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = expression("Chlorophyll b (µg/cm"^"2"*")"), fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

chlb 

# save the plot
ggsave(chlb, 
       filename = "chlb.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')


# chlorophyll c
chlc <- ggplot(data = biofilm_endpoints, 
               mapping = aes(x = factor(week), 
                             y = chlc, 
                             fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = expression("Chlorophyll c (µg/cm"^"2"*")"), fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

chlc 

# save the plot
ggsave(chlc, 
       filename = "chlc.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')

# AFDW --------------------------------------------------------------------

# mean afdw within weeks
mean_afdw_after_colonisation <- biofilm_endpoints %>%
  group_by(week, treatment) %>%
  summarize(mean_afdw_after_colonisation = mean(afdw_after_colonisation))

print(mean_afdw_after_colonisation)

# statistical differences within weeks
# week 1

biofilm_endpoints_week1 %>%
  summarize(W = shapiro.test(afdw_after_colonisation)$statistic,
            p.value = shapiro.test(afdw_after_colonisation)$p.value)

biofilm_endpoints_week1 %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw_after_colonisation, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw_after_colonisation, g = treatment)$p.value
  )


# week 2
biofilm_endpoints_week2 %>%
  summarize(W = shapiro.test(afdw_after_colonisation)$statistic,
            p.value = shapiro.test(afdw_after_colonisation)$p.value)

biofilm_endpoints_week2 %>%
  summarize(
    kruskal_statistic = kruskal.test(afdw_after_colonisation, g = treatment)$statistic,
    kruskal_p_value = kruskal.test(afdw_after_colonisation, g = treatment)$p.value
  )


# week 3
biofilm_endpoints_week3 %>%
  summarize(W = shapiro.test(afdw_after_colonisation)$statistic,
            p.value = shapiro.test(afdw_after_colonisation)$p.value)

biofilm_endpoints_week3 %>%
  summarize(
    statistic = bartlett.test(afdw_after_colonisation, treatment)$statistic,
    p.value = bartlett.test(afdw_after_colonisation, treatment)$p.value
  )

anova_week3_afdw_after_colonisation <-
  aov(
    biofilm_endpoints_week3 $ afdw_after_colonisation ~ biofilm_endpoints_week3 $ treatment
  )
summary(anova_week3_afdw_after_colonisation)

TukeyHSD(aov(anova_week3_afdw_after_colonisation))


# Plots AFDW --------------------------------------------------------------

#descriptive statistics
ggplot(biofilm_endpoints, aes(x = afdw_after_colonisation)) +
  geom_histogram(bins = 36)

# boxplot
afdw_after_colonisation_boxplot <- ggplot(data = biofilm_endpoints, 
               mapping = aes(x = factor(week), 
                             y = afdw_after_colonisation, 
                             fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = expression("Ash free dry weight (mg/cm"^"2"*")"), fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3"))

afdw_after_colonisation_boxplot 

# save the plot
ggsave(afdw_after_colonisation_boxplot, 
       filename = "afdw_after_colonisation.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')
