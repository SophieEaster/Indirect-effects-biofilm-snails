# 4 Feeding rate snails and ash free dry weight after feeding

-   NA's due to missing snail and zero values of feeding rate due to no measurable feeding in respective week

-   Unit feeding rate: µg / day and unit ash free dry weight (afdw): mg / cm²

-   Column "used_for_FA_analysis" highlights the ten snails used for FA analyses for each treatment

```{r}
# insert data
afdw_feeding_rate <-
  read.table(file.path("data\\feeding_assay_afdw_feeding_rate.txt"),
             header = TRUE, sep = "\t")

View(afdw_feeding_rate)
```

## 4.1 Feeding rate snails

#### Histogram feeding rate
```{r}
hist(afdw_feeding_rate $ feeding_rate) 
skewness(afdw_feeding_rate $ feeding_rate, na.rm = TRUE) 
shapiro.test(afdw_feeding_rate $ feeding_rate)
```
Feeding rate variable not normal distributed, right skewed and skewness > 2 -\> log10 transformation for model

Check for outliers:
```{r}
#Run function one time, after that you can use the shortcut outlierKD (your subset, variable)
outlierKD_feeding_rate <- function(dt, id_col = "treatment", start_col = 6, end_col = 6) {
  # Initialize an empty list to store results for each column
  outlier_results <- list()
  
  
  for (i in start_col:end_col) {
    col_name <- colnames(dt)[i]  # Get the column name
    
    # Extract the column values
    var_name <- dt[[col_name]]
    
    # Check if the column is numeric
    if (is.numeric(var_name)) {
      na1 <- sum(is.na(var_name))
      m1 <- mean(var_name, na.rm = TRUE)
      
      # Identify outliers using boxplot stats
      outlier <- boxplot.stats(var_name)$out  
      
      # Find rows with outliers
      outlier_rows <- which(var_name %in% outlier)
      
      # If no outliers, skip the rest of the loop
      if (length(outlier_rows) == 0) next
      
      # Get replicate values for the outliers
      id_values <- dt[[id_col]][outlier_rows]
      
      # Mean of outliers
      mo <- mean(outlier)
      
      # Collect the outlier information for this column
      outlier_results[[col_name]] <- data.frame(
        Column = col_name,
        Row = outlier_rows,
        Replicate = id_values,
        Outlier = outlier
      )
      
      cat("\nOutliers for column:", col_name, "\n")
      cat("Outliers identified:", length(outlier), "\n")
      cat("Proportion (%) of outliers:", 
          round(length(outlier) / sum(!is.na(var_name)) * 100, 1), "\n")
      cat("Mean of the outliers:", round(mo, 2), "\n")
      cat("Mean with outliers:", round(m1, 2), "\n")
    }
  }
  
  # Combine all outlier information into one data frame
  if (length(outlier_results) > 0) {
    combined_outliers <- do.call(rbind, outlier_results)
  } else {
    combined_outliers <- data.frame()  # If no outliers were found, return an empty data frame
  }
  
  # Return combined data frame with outliers information
  return(combined_outliers)
}



# Call the function and store the results
outlier_info_combined_feeding_rate <- outlierKD_feeding_rate(afdw_feeding_rate, id_col = "treatment", start_col = 6, end_col = 6)

# Print the combined outlier information
print(outlier_info_combined_feeding_rate)
```


```{r}
# Assuming outlier_info_combined is the combined outlier data frame
outlier_info_combined_feeding_rate <- outlier_info_combined_feeding_rate  # If not already assigned

# Count the number of outliers for each replicate
outlier_counts_feeding_rate <- table(outlier_info_combined_feeding_rate$Replicate)

# Convert to a data frame for better readability
outlier_counts_feeding_rate_df <- as.data.frame(outlier_counts_feeding_rate)
colnames(outlier_counts_feeding_rate_df) <- c("Replicate", "Outlier_Count")

View(outlier_counts_feeding_rate_df)
```

Remove outliers:
```{r}
# Extract the row numbers of the outliers
outlier_rows_feeding_rate <- outlier_info_combined_feeding_rate$Row

# Remove the outlier rows from the dataset
afdw_feeding_rate_cleaned_1 <- afdw_feeding_rate[-outlier_rows_feeding_rate, ]

```
Check histogram again:
```{r}
hist(afdw_feeding_rate_cleaned_1 $ feeding_rate) 
skewness(afdw_feeding_rate_cleaned_1 $ feeding_rate, na.rm = TRUE) 
shapiro.test(afdw_feeding_rate_cleaned_1 $ feeding_rate)
```
Still not normal distributed and right skewed, but skewness now 1.3, therefore no data transformation for later model.

#### Calculate mean feeding rate for each week

```{r}
mean_feeding_rate <- afdw_feeding_rate_cleaned_1 %>%
  group_by(week, treatment) %>%
  summarize(mean_feeding_rate = mean(feeding_rate, na.rm = TRUE))

print(mean_feeding_rate)
```

Now calculate mean feeding for each treatment:

```{r}
mean_feeding_rate_treatment <- afdw_feeding_rate_cleaned_1 %>%
  group_by(treatment) %>%
  summarize(mean_feeding_rate = mean(feeding_rate, na.rm = TRUE))

print(mean_feeding_rate_treatment)
```

#### Model feeding rate

```{r}
# week, treatment and beaker as factor
afdw_feeding_rate_cleaned_1 $ week <- factor(afdw_feeding_rate_cleaned_1 $ week)
afdw_feeding_rate_cleaned_1 $ treatment <- factor(afdw_feeding_rate_cleaned_1 $ treatment)
afdw_feeding_rate_cleaned_1 $ beaker <- factor(afdw_feeding_rate_cleaned_1 $ beaker)

# response variable: feeding_rate
# fixed effects: week + treatment + week:treatment
# random effect with constant intercept: (1|beaker)
# first check for NAs
any(is.na(afdw_feeding_rate_cleaned_1$feeding_rate))

feeding_rate_lmer <- lmer(feeding_rate ~ week + treatment + week:treatment + (1|beaker), data = afdw_feeding_rate_cleaned_1)

par(mfrow = c(2, 2))

plot(feeding_rate_lmer) # point cloud looks evenly distributed --> normal distribution

```

```{r}
qreference(residuals(feeding_rate_lmer)) # normality ok

drop1(feeding_rate_lmer, test = 'Chisq')
# interaction of timepoint (week) and treatment highly significant
# now check at which timepoints is a difference between treatments
anova(feeding_rate_lmer)

# conduct Tukey HSD-test for pairwise comparisons
feeding_rate_tukey <- lsmeans(feeding_rate_lmer, ~ treatment | week)
feeding_rate_tukey

feeding_rate_tukey1 <- contrast(feeding_rate_tukey, 'pairwise')
feeding_rate_tukey1

feeding_rate_summary <- summary(feeding_rate_tukey1, infer = TRUE, adjust = 'mvt')
feeding_rate_summary

```

Differences: 
week 1: Cip-Pro and Mix_pro
week 2: Mix-NC and NC-Pro
week 3: none

```{r}
feeding_rate_inter1 <- lsmeans(feeding_rate_lmer, ~ week |  treatment)

feeding_rate_inter_dunn1 <- contrast(feeding_rate_inter1, 'pairwise')

(feeding_rate_inter_dunn_sum1 <- summary(feeding_rate_inter_dunn1, infer = TRUE, adjust = 'mvt'))
```

Weeks for all treatments similar, only difference for treatment NC: week 2 & 3 and for treatment Pro: week 1 & 3

#### Boxplot feeding rate after 7/14/21 days

```{r}
feeding_rate_boxplot <- ggplot(data = afdw_feeding_rate_cleaned_1, 
       mapping = aes(x = factor(week), 
                     y = feeding_rate, 
                     fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3")) +
  labs(y = expression("Feeding rate (µg/day)"))

feeding_rate_boxplot

# save the plot
ggsave(feeding_rate_boxplot, 
       filename = "feeding_rate_boxplot.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\Bottom_up_effects_biofilm_grazer\\output')


```

## 4.2 AFDW after feeding

Calculate mean of afdw after feeding for each week

```{r}
mean_afdw_after_feeding <- afdw_feeding_rate %>%
  group_by(week, treatment) %>%
  summarize(mean_afdw_after_feeding = mean(afdw))

print(mean_afdw_after_feeding)
```
Calculate mean of afdw after feeding for each treatment
```{r}
mean_afdw_after_feeding_treatment <- afdw_feeding_rate %>%
  group_by(treatment) %>%
  summarize(mean_afdw_after_feeding = mean(afdw))

print(mean_afdw_after_feeding_treatment)
```


#### Histogram afdw after feeding

```{r}
hist(afdw_feeding_rate $ afdw) 
skewness(afdw_feeding_rate $ afdw, na.rm = TRUE) 
shapiro.test(afdw_feeding_rate $ afdw)
```

AFDW after feeding not normal distributed, skewness = 1.05, no data transformation.

#### Model afdw after feeding

```{r}
# response variable: afdw
# fixed effects: week + treatment + week:treatment
# random effect with constant intercept: (1|beaker)
# remove NAs
afdw_feeding_rate <- na.omit(afdw_feeding_rate)

afdw_after_feeding_lmer <- lmer(afdw ~ week + treatment + week:treatment + (1|beaker), data = afdw_feeding_rate)

#----- new try without random effect
model_afdw2 <- lm(afdw ~ week + treatment + week:treatment, data = afdw_feeding_rate)

drop1(model_afdw2, test = "Chisq")

anova(model_afdw2)
#----------
plot(afdw_after_feeding_lmer) # point cloud looks evenly distributed --> normal distribution

qreference(residuals(afdw_after_feeding_lmer)) # normality

drop1(afdw_after_feeding_lmer, test = 'Chisq')
# interaction of timepoint (week) and treatment highly significant
# now check at which timepoints is a difference between treatments
anova(afdw_after_feeding_lmer)

# conduct Tukey HSD-test for pairwise comparisons
afdw_tukey <- lsmeans(afdw_after_feeding_lmer, ~ treatment | week)
afdw_tukey

afdw_tukey1 <- contrast(afdw_tukey, 'pairwise')
afdw_tukey

afdw_summary <- summary(afdw_tukey1, infer = TRUE, adjust = 'mvt')
afdw_summary
```

Here treatments quite different within weeks. Compare weeks:

```{r}
afdw_after_feeding_inter1 <- lsmeans(afdw_after_feeding_lmer, ~ week |  treatment)

afdw_after_feeding_dunn1 <- contrast(afdw_after_feeding_inter1, 'pairwise')

(afdw_after_feeding_dunn_sum1 <- summary(afdw_after_feeding_dunn1, infer = TRUE, adjust = 'mvt'))
```

Week 1 & 2 for treatments Mix & Pro similar. Everything else shows differences.

#### Boxplot afdw after feeding of 7/14/21 days

```{r}
afdw_after_feeding_boxplot <- ggplot(data = afdw_feeding_rate, 
       mapping = aes(x = factor(week), 
                     y = afdw, 
                     fill = factor(treatment, levels = c("NC", "Cip", "Pro", "Mix")))) +
  geom_boxplot() +
  scale_fill_manual(values = c("NC" = "#009E73", "Cip" = "#56B4E9", "Pro" = "#E69F00", "Mix" = "#F0E442")) + 
  scale_color_okabe_ito() + 
  theme_bw() +
  labs(x = NULL , fill = "Treatment") +
  scale_x_discrete(labels = c("Week 1", "Week 2", "Week 3")) +
  labs(y = expression("AFDW (mg/cm"^"2"*")"))

afdw_after_feeding_boxplot

# save the plot
ggsave(afdw_after_feeding_boxplot, 
       filename = "afdw_after_feeding_boxplot.png", 
       path = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\Bottom_up_effects_biofilm_grazer\\output')
```