# AFDW afer 7 days of feeding
library('ggplot2')
library(ggsignif)
library(car)
library(data.table)
library(dplyr)
library(rstatix)
# Insert data
AFDW_after_feeding <- read.table(file = 'C:\\Users\\sophi\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Snail experiment\\AFDW_after_feeding.txt', 
                               header = TRUE, sep = "\t")


AFDW_after_feeding <- na.omit(AFDW_after_feeding)
View(AFDW_after_feeding)

# One-way ANOVA
# AFDW as a function of the treatments, subsetted into the 3 weeks
# Subset into weeks and statistics
Week1 <- subset(AFDW_after_feeding, Week == "Week1", select = c(Sample, Treatment, AFDW))
View(Week1)
# Median Week 1
median_1 <- Week1 %>%
  group_by(Treatment) %>%
  summarize(Week1 = median(AFDW))
print(median_1)
# Significance?
shapiro.test(Week1$AFDW)
kruskal.test(Week1$AFDW ~ Week1$Treatment)
# p=4.792e-06
# post hoc test
dunn_test(Week1, AFDW ~ Treatment, p.adjust.method = "holm")
# NC-Cip: p=0.000944***
# NC-Pro: p=0.000125***

#Week2
Week2 <- subset(AFDW_after_feeding, Week == "Week2", select = c(Sample, Treatment, AFDW))
View(Week2)
# Median Week 2
median_2 <- Week2 %>%
  group_by(Treatment) %>%
  summarize(Week2 = median(AFDW))
print(median_2)
# Significance?
shapiro.test(Week2$AFDW)
kruskal.test(Week2$AFDW ~ Week2$Treatment)
# p=1.014e-09
# post hoc test
dunn_test(Week2, AFDW ~ Treatment, p.adjust.method = "holm")
# NC-Cip: p=ns
# NC-Mix: p=0.000000631****
# NC-Pro: p=0.0198*

#Week3
Week3 <- subset(AFDW_after_feeding, Week == "Week3", select = c(Sample, Treatment, AFDW))
View(Week3)
# Median Week 3
median_3 <- Week3 %>%
  group_by(Treatment) %>%
  summarize(Week3 = median(AFDW))
print(median_3)
# Significance?
shapiro.test(Week3$AFDW)
bartlett.test(AFDW ~ Treatment, Week3)
anova_week3 <- aov(Week3$AFDW~Week3$Treatment)
summary(anova_week3)
# p=2.26e-11 ***
# post hoc test
TukeyHSD(anova_week3)
# NC-Mix: p=0.0000144
# NC-Pro: p=0.02121

#--------------------------------------------------------------------
# check:  https://www.scribbr.com/statistics/anova-in-r/
# Two-way ANOVA
# AFDW as a function of the treatment and the Week
two.way <- aov(AFDW ~ Treatment + Week, data=AFDW_after_feeding)
summary(two.way)

#--------------------------------------------------------------------

# Check for homoscedasticity

par(mfrow=c(2,2))
plot(two.way)
par(mfrow=c(1,1))

#post hoc text
tukey.two.way <- TukeyHSD(two.way)
tukey.two.way


#-------------------------------------------------------------------
# Plot results
tukey.plot.aov<-aov(AFDW ~ Treatment:Week, data=AFDW_after_feeding)

tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)

new_dataframe <- AFDW_after_feeding %>%
  group_by(Week, Treatment) %>%
  summarise(
    AFDW = median(AFDW)
  )
# Next, add the group labels as a new variable in the data frame.
#### to be continued
