# graphs
# packages
library(stats)
library(rstatix)
library(dplyr)
library(tidyverse)
library(tidyr)
library(broom)
library(dunn.test)

# insert data 
# 1) feeding rate
feeding_rate <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\feeding_rate.txt',
             header = TRUE, sep = "\t")

feeding_rate <- na.omit(feeding_rate)

View(feeding_rate)

# 2) growth snails
growth_snails <- 
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\growth_snails.csv',
             header = TRUE, sep = ",")

filtered_growth <- growth_snails %>%
  filter(!treatment %in% c("Culture"))


# combine tables, row sizes must be compatible when column-binding
feeding_and_growth <-
  bind_cols(feeding_rate, growth_snails)
