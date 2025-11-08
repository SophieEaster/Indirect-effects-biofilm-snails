# Dryweight snails
library('ggplot2')
library(ggsignif)
library(car)
library(data.table)
library(dplyr)
# Insert data
Dryweight_snails <- read.table(file = 'C:\\Users\\sophi\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Snail experiment\\Dryweight_snails.txt', 
                             header = TRUE, sep = "\t")


View(Dryweight_snails)
# Median
median_dw_snails <- Dryweight_snails %>%
  group_by(Treatment) %>%
  summarize(median_dw_snails = median(dw_snail))
print(median_dw_snails)
# Significance?
shapiro.test(Dryweight_snails$dw_snail)
bartlett.test(dw_snail ~ Treatment, Dryweight_snails)
anova_dw_snail <- aov(Dryweight_snails$dw_snail~Dryweight_snails$Treatment)
summary(anova_dw_snail)
# p=0.000571***
#post hoc test
TukeyHSD(aov(Dryweight_snails$dw_snail~Dryweight_snails$Treatment))
# NC-Cip: p=0.0351
# NC-Mix: p=0.00797 
# NC-Culture: p=0.0002354 
# NC-Pro: p=0.0873