# Dryweight of snails

# packages
library(base)
library(stats)
library(ggplot2)
# data
dryweight_snails <- read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\dryweight_snails.txt',
                               header = TRUE, sep = "\t")

View(dryweight_snails)

# mean of each treatment
mean_dryweight_snails <- dryweight_snails %>%
              group_by(treatment) %>%
              summarize(mean_dryweight_snails = mean(dryweight))

print(mean_dryweight_snails)

shapiro.test(dryweight_snails$dryweight)
bartlett.test(dryweight ~ treatment, dryweight_snails)
anova_dryweight_snails <- aov(dryweight_snails$dryweight ~ dryweight_snails$treatment)
summary(anova_dryweight_snails)
TukeyHSD(aov(anova_dryweight_snails))

# calculate growth
# mean baseline dryweight of 10 lyophilized culture snails = 0.452
growth_snails <- dryweight_snails %>%
  dplyr::mutate(growth_ng_dw_per_day = ((dryweight - 0.452) / 7) * 1000)

# save for usage in graphs.R
write.csv(growth_snails, file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\growth_snails.csv'
)

# statistics
filtered_growth <- growth_snails %>%
  filter(!treatment %in% c("Culture"))

filtered_growth %>%
  get_summary_stats(dryweight, growth_ng_dw_per_day)


ggplot(filtered_growth, aes(x = dryweight)) +
               geom_histogram(binwidth = 0.05)
               

mean_growth_snails <- filtered_growth %>%
  group_by(treatment) %>%
  summarize(mean_growth_snails = mean(growth_ng_dw_per_day))
  
print(mean_growth_snails)  

shapiro_growth <- filtered_growth %>%
  summarize(W = shapiro.test(growth_ng_dw_per_day)$statistic,
            p.value = shapiro.test(growth_ng_dw_per_day)$p.value)

print(shapiro_growth)

bartlett_growth <- filtered_growth %>%
  summarize(
    statistic = bartlett.test(growth_ng_dw_per_day, treatment)$statistic,
    p.value = bartlett.test(growth_ng_dw_per_day, treatment)$p.value
  )

print(bartlett_growth)

anova_growth_snails <- aov(growth_ng_dw_per_day ~ treatment, data = filtered_growth)
summary(anova_growth_snails)

TukeyHSD(aov(anova_growth_snails))
