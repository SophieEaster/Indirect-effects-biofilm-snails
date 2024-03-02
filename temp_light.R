# Temperature and light 

# Feeding assay -----------------------------------------------------------

# 1) Feeding assay

# packages
library(stats)
library(rstatix)
library(dplyr)
library(tidyverse)
library(tidyr)
library(broom)
library(dunn.test)
library(patchwork)
library(ggplot2)
library(ggokabeito)

# insert data
feeding_assay_temp_light <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\feeding_assay_temp_light.txt',
             header = TRUE, sep = "\t")

View(feeding_assay_temp_light)

# average temperature
mean(feeding_assay_temp_light $ temperature)
# 13.27927 °C

# mean light intensity
# filter only daytime
average_light_intensity_feeding <- feeding_assay_temp_light %>%
  dplyr::filter(light_intensity > 0) %>%
  summarize(average_light_intensity_feeding = mean(light_intensity))

print(average_light_intensity_feeding)
# 330.3552 lux

# temperature plot
# add a small amount of random noise to the y-values
feeding_assay_temp_light $ temperature <- jitter(feeding_assay_temp_light $ temperature)

t1 <- ggplot(
  data = feeding_assay_temp_light,
  mapping = aes( x = day, y = temperature, color = time_of_day)) +
  geom_point(alpha = 0.5) +  # Individual measurements with transparency
  geom_smooth() +
  labs(
    title = "Temperature profile during 21-day feeding assay",
    x = "Day", y = "Temperature (°C)"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = 0:21) +
  scale_color_continuous(
    name = "Time of Day",  
    breaks = c(0, 5, 10, 15, 20),
    labels = c("0:00 am", "05:00 am", "10:00 am", "03:00 pm", "08:00 pm"),
    limits = c(0, 22),
  )

t1 

# save the plot
ggsave(t1, 
       filename = "temperature_feeding_assay.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')


# light intensity plot
# Filter the data to include only rows where light_intensity > 0
filtered_data1 <- feeding_assay_temp_light %>% 
  filter(light_intensity > 0)

# add a small amount of random noise to the y-values
filtered_data1 $ light_intensity <- jitter(filtered_data1 $ light_intensity)

l1 <- ggplot(
  data = filtered_data1,
  mapping = aes( x = day, y = light_intensity)) +
  geom_point(size = 1) +  # Individual measurements with transparency
  geom_smooth(color = "#E69F00", linewidth = 1) +
  labs(
    title = "Light intensity profile during 21-day feeding assay",
    x = "Day", y = "Light intensity (lux)"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = 0:21) 


l1 

# save the plot
ggsave(l1, 
       filename = "light_intensity_feeding_assay.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')


# Biofilm colonisation ----------------------------------------------------
# insert data
colonisation_temp_light <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\colonisation_temp_light.txt',
             header = TRUE, sep = "\t")

View(colonisation_temp_light)

# average temperature
mean(colonisation_temp_light $ temperature_colonisation)
# 18.28751 °C

# mean light intensity
# filter only daytime
average_light_intensity_colonisation <- colonisation_temp_light %>%
  dplyr::filter(light_intensity_colonisation > 0) %>%
  summarize(average_light_intensity_colonisation = mean(light_intensity_colonisation))

print(average_light_intensity_colonisation)
# 14019.26 lux

# temperature plot
# add a small amount of random noise to the y-values
colonisation_temp_light $ temperature_colonisation <- jitter(colonisation_temp_light $ temperature_colonisation)

t2 <- ggplot(
  data = colonisation_temp_light,
  mapping = aes( x = day, y = temperature_colonisation, color = time_of_day)) +
  geom_point(alpha = 0.5) +  # Individual measurements with transparency
  geom_smooth() +
  labs(
    title = "Temperature profile during biofilm colonisation",
    x = "Day", y = "Temperature (°C)"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 77, by = 5)) +
  scale_color_continuous(
    name = "Time of Day",  
    breaks = c(0, 5, 10, 15, 20),
    labels = c("0:00 am", "05:00 am", "10:00 am", "03:00 pm", "08:00 pm"),
    limits = c(0, 22),
  )

t2 

# save the plot
ggsave(t2, 
       filename = "temperature_colonisation.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')

# light intensity plot
# Filter the data to include only rows where light_intensity > 0
filtered_data2 <- colonisation_temp_light %>% 
  filter(light_intensity_colonisation > 0)

# add a small amount of random noise to the y-values
filtered_data2 $ light_intensity_colonisation <- jitter(filtered_data2 $ light_intensity_colonisation)

l2 <- ggplot(
  data = filtered_data2,
  mapping = aes( x = day, y = light_intensity_colonisation)) +
  geom_point(size = 1) +  # Individual measurements with transparency
  geom_smooth(color = "#E69F00", linewidth = 1) +
  labs(
    title = "Light intensity profile during biofilm colonisation",
    x = "Day", y = "Light intensity (lux)"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 77, by = 5)) 


l2 

# save the plot
ggsave(l2, 
       filename = "light_intensity_colonisation.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\output')

