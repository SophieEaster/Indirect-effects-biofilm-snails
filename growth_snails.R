# Dryweight snails / growth of snails / combination of snail growth and feeding rate and PUFA content

# packages
library(stats)
library(rstatix)
library(dplyr)
library(tidyverse)
library(stringr)
library(tidyr)
library(broom)
library(dunn.test)
library(ggplot2)
library(ggokabeito)
library(ggpubr)

# insert data
dryweight_snails <- read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\dryweight_snails.txt',
                               header = TRUE, sep = "\t")

View(dryweight_snails)


# Dryweight snails -------------------------------------------------

# mean dryweight of each treatment, unit = mg
mean_dryweight_snails <- dryweight_snails %>%
              group_by(treatment) %>%
              summarize(mean_dryweight_snails = mean(dryweight_mg))

print(mean_dryweight_snails)

shapiro.test(dryweight_snails $ dryweight_mg)
bartlett.test(dryweight_mg ~ treatment, dryweight_snails)
anova_dryweight_snails <- aov(dryweight_snails $ dryweight_mg ~ dryweight_snails $ treatment)
summary(anova_dryweight_snails)
TukeyHSD(aov(anova_dryweight_snails))


# Growth snails -----------------------------------------------------------
# calculate growth
# mean baseline dryweight of 10 lyophilized culture snails = 0.452 mg
growth_snails <- dryweight_snails %>%
  dplyr::mutate(growth_µg_dw_per_day = ((dryweight_mg - 0.452) / 7) * 1000) # multiply by 1000 to obtain a growth unit of µg dryweight per day

# statistics
filtered_growth <- growth_snails %>%
  filter(!treatment %in% c("Culture"))

filtered_growth %>%
  get_summary_stats(dryweight_mg, growth_µg_dw_per_day)


ggplot(filtered_growth, aes(x = growth_µg_dw_per_day)) +
               geom_histogram(binwidth = 10)
               

mean_growth_snails <- filtered_growth %>%
  group_by(treatment) %>%
  summarize(mean_growth_snails = mean(growth_µg_dw_per_day))
  
print(mean_growth_snails)  

shapiro_growth <- filtered_growth %>%
  summarize(W = shapiro.test(growth_µg_dw_per_day)$statistic,
            p.value = shapiro.test(growth_µg_dw_per_day)$p.value)

print(shapiro_growth)

bartlett_growth <- filtered_growth %>%
  summarize(
    statistic = bartlett.test(growth_µg_dw_per_day, treatment)$statistic,
    p.value = bartlett.test(growth_µg_dw_per_day, treatment)$p.value
  )

print(bartlett_growth)

anova_growth_snails <- aov(growth_µg_dw_per_day ~ treatment, data = filtered_growth)
summary(anova_growth_snails)

TukeyHSD(aov(anova_growth_snails))


# Plot growth snails ------------------------------------------------------
# box plot growth after 21-day feeding assay


growth_snails_boxplot <- ggplot(data = filtered_growth, 
                                mapping = aes(x = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")), 
                                              y = growth_µg_dw_per_day, 
                                              fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442"),
                    breaks = c("NC", "Cip", "Pro", "Mix")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , y = "Growth (µg dryweight / day)", fill = "Treatment") +
  scale_x_discrete(labels = c("NC", "Cip", "Pro", "Mix")) +
  labs(title = expression(paste("Growth of ", italic("Potamopyrgus antipodarum"), " after 21-day feeding assay"))) +
  annotate("text", x = 1.5, y = max(filtered_growth$growth_µg_dw_per_day) + 3, 
           label = "*", size = 6) +
  geom_segment(aes(x = 1, xend = 2, y = max(growth_µg_dw_per_day) + 2, yend = max(growth_µg_dw_per_day) + 2),
               lineend = "round", linewidth = 0.8) +
  annotate("text", x = 2.5, y = max(filtered_growth$growth_µg_dw_per_day) + 6, 
           label = "**", size = 6) +
  geom_segment(aes(x = 1, xend = 4, y = max(growth_µg_dw_per_day) + 5, yend = max(growth_µg_dw_per_day) + 5),
               lineend = "round", linewidth = 0.8) 


growth_snails_boxplot

# save the plot
ggsave(growth_snails_boxplot, 
       filename = "growth_snails_boxplot.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')


# Combination growth & feeding rate ---------------------------------------

# insert feeding rate data set
afdw_feeding_rate <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\feeding_assay_afdw_feeding_rate.txt',
             header = TRUE, sep = "\t")

View(afdw_feeding_rate)

# calculate the total feeding rate after 21 days for the 40 snails used in the feeding assay, afdw not needed here
# first create new data frame only with snails used in the feeding assay
total_feeding_rate <-
  afdw_feeding_rate %>%
  filter(startsWith(used_in_feeding_assay, "yes")) %>%
  select(-afdw)

# now add the feeding rates of the 3 weeks
total_feeding_per_snail <- total_feeding_rate %>%
  group_by(used_in_feeding_assay) %>%
  summarize(total_feeding = sum(feeding_rate))

print(total_feeding_per_snail)

# Remove the prefix "yes_"
total_feeding_per_snail <- total_feeding_per_snail %>%
  mutate(used_in_feeding_assay = str_remove(used_in_feeding_assay, "^yes_"))

print(total_feeding_per_snail)

# Reorder data frame
# Define new order of rows
new_order <- c("NC1", "NC2", "NC3", "NC4", "NC5", "NC6", "NC7", "NC8", "NC9", "NC10",
               "Pro1", "Pro2", "Pro3", "Pro4", "Pro5", "Pro6", "Pro7", "Pro8", "Pro9", "Pro10",
               "Cip1", "Cip2", "Cip3", "Cip4", "Cip5", "Cip6", "Cip7", "Cip8", "Cip9", "Cip10",
               "Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6", "Mix7", "Mix8", "Mix9", "Mix10")

# Reorder the rows based on the 'used_in_feeding_assay' column
total_feeding_per_snail <- total_feeding_per_snail[match(new_order, total_feeding_per_snail $ used_in_feeding_assay), ]

print(total_feeding_per_snail)

# Join tables "filtered_growth" and "total_feeding_per_snail"
feeding_vs_growth <- bind_cols(filtered_growth, total_feeding_per_snail)

### # Try out model

# statistical analysis

hist(feeding_vs_growth$growth_µg_dw_per_day) # normal distributed

qqnorm(feeding_vs_growth$growth_µg_dw_per_day, datax=TRUE)
qqline(feeding_vs_growth$growth_µg_dw_per_day, datax=TRUE)

# model building
# linear mixed model with 


# Assuming you want to perform the linear regression only on the "NC" treatment
treatment_subset <- "NC"  # Specify the treatment you want to include

# Subset the data to include only rows with the specified treatment
subset_data <- subset(feeding_vs_growth, treatment == treatment_subset)

# Linear regression
# explanatory variables: feeding and treatment, response variable: growth
lm_feeding_vs_growth = lm(growth_µg_dw_per_day ~ total_feeding, data = subset_data)
summary(lm_feeding_vs_growth)
