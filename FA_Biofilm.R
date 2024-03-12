# Fatty acid analyses for biofilm 
# load packages
library(vegan)
library(data.table)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(patchwork)
library(dplyr)
library(stringr)
library(tidyr)
library(shades)
library(ggpubr)
library(base)
library(ggokabeito)
library(gridExtra)
library(multcomp)
library(car)
library(data.table)
library(asbio)
library(plotrix)
library(broom)
library(lmtest)
library(boot)

# FA_Biofilm_absolute -----------------------------------------------------

# insert data
FA_Biofilm_absolute <-
  read.table(file = 'C:\\Users\\Oster Test\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Indirect-effects-biofilm-snails\\data\\FA_Biofilm_absolute.txt',
             header = TRUE, sep = "\t")


View(FA_Biofilm_absolute)



# NMDS biofilm absolute ---------------------------------------------------
# add prefix (week) 1,2 or 3 to treatments 
FA_Biofilm_absolute <- FA_Biofilm_absolute %>%
  mutate(treatment = if_else(week == "1", paste("1_", treatment, sep = ""), treatment))

FA_Biofilm_absolute <- FA_Biofilm_absolute %>%
  mutate(treatment = if_else(week == "2", paste("2_", treatment, sep = ""), treatment))

FA_Biofilm_absolute <- FA_Biofilm_absolute %>%
  mutate(treatment = if_else(week == "3", paste("3_", treatment, sep = ""), treatment))


#change to wide format
FA_Biofilm_absolute_wide <-
  pivot_wider(FA_Biofilm_absolute,
              names_from = fa,
              values_from = fa_conc)

# square root transformation
FA_Biofilm_absolute_wide <- sqrt(FA_Biofilm_absolute_wide[, -c(1, 2, 3)])
View(FA_Biofilm_absolute_wide)

#set NA to 0 
FA_Biofilm_absolute_wide[is.na(FA_Biofilm_absolute_wide)] <- 0
View(FA_Biofilm_absolute_wide)

#meta MDS
NMDS_Biofilm_absolute <-
  metaMDS(
    FA_Biofilm_absolute_wide,
    autotransform = FALSE,
    k = 2 ,
    distance = "bray",
    na.rm = TRUE
  )

NMDS_Biofilm_absolute$stress

#add treatment column
colnames(FA_Biofilm_absolute_wide) <-
  c(
    "C08:0",
    "C10:0",
    "C11:0",
    "C12:0",
    "C13:0",
    "C14:0",
    "C14:1",
    "C15:0",
    "C15:1",
    "C16:0",
    "C16:1",
    "C17:0",
    "C17:1",
    "C18:0",
    "C18:1",
    "C18:2(n-6t)",
    "C18:2(n-6c)",
    "C18:3(n-6)",
    "C18:3(n-3)",
    "C18:4(n-3)",
    "C20:0",
    "C20:1",
    "C20:2(n-6)",
    "C20:3(n-6)",
    "C20:4(n-6)",
    "C20:5(n-3)",
    "C21:0",
    "C22:1",
    "C22:0",
    "C22:5(n-3)",
    "C22:6(n-3)",
    "C24:0",
    "C24:1"
  )

#environmental factors and vectors
# "ba" = biofilm absolute
FA_variables.SAFAs_ba <-
  envfit(
    NMDS_Biofilm_absolute,
    subset(
      FA_Biofilm_absolute_wide,
      select = c(
        "C08:0",
        "C10:0",
        "C11:0",
        "C12:0",
        "C13:0",
        "C14:0",
        "C15:0",
        "C16:0",
        "C17:0",
        "C18:0",
        "C20:0",
        "C21:0",
        "C22:0",
        "C24:0"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.MUFAs_ba <-
  envfit(
    NMDS_Biofilm_absolute,
    subset(
      FA_Biofilm_absolute_wide,
      select = c(
        "C14:1",
        "C15:1",
        "C16:1",
        "C17:1",
        "C18:1",
        "C20:1",
        "C22:1",
        "C24:1"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.PUFAs_ba <-
  envfit(
    NMDS_Biofilm_absolute,
    subset(
      FA_Biofilm_absolute_wide,
      select = c(
        "C18:2(n-6t)",
        "C18:2(n-6c)",
        "C18:3(n-6)",
        "C18:3(n-3)",
        "C18:4(n-3)",
        "C20:2(n-6)",
        "C20:3(n-6)",
        "C20:4(n-6)",
        "C20:5(n-3)",
        "C22:5(n-3)",
        "C22:6(n-3)"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

#data frame
data.scores_biofilm_absolute = as.data.frame(scores(NMDS_Biofilm_absolute)$sites)

Treatments <-
  FA_Biofilm_absolute$treatment[FA_Biofilm_absolute$fa == "C16:0"]

data.scores_biofilm_absolute$treatment = Treatments

# create weeks as subgroups for plotting

data.scores_biofilm_absolute <- data.scores_biofilm_absolute %>%
  mutate(subgroup = case_when(
    str_detect(treatment, "^1_") ~ "Week1",
    str_detect(treatment, "^2_") ~ "Week2",
    str_detect(treatment, "^3_") ~ "Week3",
    TRUE ~ "Other"
  ))

# plot NMDS

absol_Biofilm = ggplot(data = data.scores_biofilm_absolute, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    data = data.scores_biofilm_absolute,
    aes(colour = subgroup),
    size = 5,
    alpha = 0.5
  ) +
  geom_text_repel(
    aes(label = str_sub(treatment, start = 3)),
    box.padding = 0.5, # increase padding around text
    force = 10, # increase the repulsion strength
    segment.color = NA, # to remove connecting line
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    #axis.ticks = element_blank(),
    #axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 10, colour = "grey30"),
    plot.margin = margin(
      t = 0.2,
      b = 0.2,
      l = 0.2,
      r = 0.2,
      unit = "cm"
    ),
    panel.grid.major = element_line(colour = "grey95", size = 0.4),
    panel.grid.minor = element_line(colour = "grey95", size = 0.2)
  ) +
  annotate(
    "text",
    x = 0.3,
    y = -0.15,
    label = paste("stress:", NMDS_Biofilm_absolute$stress %>% round(2)),
    size = 3
  )

absol_Biofilm




# FA_Biofilm_relative -----------------------------------------------------

## hier noch umrechnen in %





# NMDS biofilm relative ---------------------------------------------------
# add prefix (week) 1,2 or 3 to treatments 
FA_Biofilm_relative <- FA_Biofilm_relative %>%
  mutate(treatment = if_else(week == "1", paste("1_", treatment, sep = ""), treatment))

FA_Biofilm_relative <- FA_Biofilm_relative %>%
  mutate(treatment = if_else(week == "2", paste("2_", treatment, sep = ""), treatment))

FA_Biofilm_relative <- FA_Biofilm_relative %>%
  mutate(treatment = if_else(week == "3", paste("3_", treatment, sep = ""), treatment))


#change to wide format
FA_Biofilm_relative_wide <-
  pivot_wider(FA_Biofilm_relative,
              names_from = fa,
              values_from = fa_rel)

# square root transformation
FA_Biofilm_relative_wide <- sqrt(FA_Biofilm_relative_wide[, -c(1, 2, 3)])
View(FA_Biofilm_relative_wide)

#set NA to 0 
FA_Biofilm_relative_wide[is.na(FA_Biofilm_relative_wide)] <- 0
View(FA_Biofilm_relative_wide)

#meta MDS
NMDS_Biofilm_relative <-
  metaMDS(
    FA_Biofilm_relative_wide,
    autotransform = FALSE,
    k = 2 ,
    distance = "bray",
    na.rm = TRUE
  )

NMDS_Biofilm_relative$stress

#add treatment column
colnames(FA_Biofilm_relative_wide) <-
  c(
    "C08:0",
    "C10:0",
    "C11:0",
    "C12:0",
    "C13:0",
    "C14:0",
    "C14:1",
    "C15:0",
    "C15:1",
    "C16:0",
    "C16:1",
    "C17:0",
    "C17:1",
    "C18:0",
    "C18:1",
    "C18:2(n-6t)",
    "C18:2(n-6c)",
    "C18:3(n-6)",
    "C18:3(n-3)",
    "C18:4(n-3)",
    "C20:0",
    "C20:1",
    "C20:2(n-6)",
    "C20:3(n-6)",
    "C20:4(n-6)",
    "C20:5(n-3)",
    "C21:0",
    "C22:1",
    "C22:0",
    "C22:5(n-3)",
    "C22:6(n-3)",
    "C24:0",
    "C24:1"
  )

#environmental factors and vectors
# "br" = biofilm relative
FA_variables.SAFAs_br <-
  envfit(
    NMDS_Biofilm_relative,
    subset(
      FA_Biofilm_relative_wide,
      select = c(
        "C08:0",
        "C10:0",
        "C11:0",
        "C12:0",
        "C13:0",
        "C14:0",
        "C15:0",
        "C16:0",
        "C17:0",
        "C18:0",
        "C20:0",
        "C21:0",
        "C22:0",
        "C24:0"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.MUFAs_br <-
  envfit(
    NMDS_Biofilm_relative,
    subset(
      FA_Biofilm_relative_wide,
      select = c(
        "C14:1",
        "C15:1",
        "C16:1",
        "C17:1",
        "C18:1",
        "C20:1",
        "C22:1",
        "C24:1"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

FA_variables.PUFAs_br <-
  envfit(
    NMDS_Biofilm_relative,
    subset(
      FA_Biofilm_relative_wide,
      select = c(
        "C18:2(n-6t)",
        "C18:2(n-6c)",
        "C18:3(n-6)",
        "C18:3(n-3)",
        "C18:4(n-3)",
        "C20:2(n-6)",
        "C20:3(n-6)",
        "C20:4(n-6)",
        "C20:5(n-3)",
        "C22:5(n-3)",
        "C22:6(n-3)"
      )
    )
    ,
    permutations = 999,
    na.rm = TRUE
  )

#data frame
data.scores_biofilm_relative = as.data.frame(scores(NMDS_Biofilm_relative)$sites)

Treatments <-
  FA_Biofilm_relative$treatment[FA_Biofilm_relative$fa == "C16:0"]

data.scores_biofilm_relative$treatment = Treatments


# create weeks as subgroups for plotting

data.scores_biofilm_relative <- data.scores_biofilm_relative %>%
  mutate(subgroup = case_when(
    str_detect(treatment, "^1_") ~ "Week1",
    str_detect(treatment, "^2_") ~ "Week2",
    str_detect(treatment, "^3_") ~ "Week3",
    TRUE ~ "Other"
  ))

# plot NMDS

rel_Biofilm = ggplot(data = data.scores_biofilm_relative, aes(x = NMDS1, y = NMDS2)) +
  geom_point(
    data = data.scores_biofilm_relative,
    aes(colour = subgroup),
    size = 5,
    alpha = 0.5
  ) +
  geom_text_repel(
    aes(label = str_sub(treatment, start = 3)),
    box.padding = 0.5, # increase padding around text
    force = 10, # increase the repulsion strength
    segment.color = NA, # to remove connecting line
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey30"),
    #axis.ticks = element_blank(),
    #axis.text = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 11,
      face = "bold",
      colour = "grey30"
    ),
    legend.text = element_text(size = 10, colour = "grey30"),
    plot.margin = margin(
      t = 0.2,
      b = 0.2,
      l = 0.2,
      r = 0.2,
      unit = "cm"
    ),
    panel.grid.major = element_line(colour = "grey95", size = 0.4),
    panel.grid.minor = element_line(colour = "grey95", size = 0.2)
  ) +
  annotate(
    "text",
    x = -0.1,
    y = -0.15,
    label = paste("stress:", NMDS_Biofilm_absolute$stress %>% round(2)),
    size = 3
  )

rel_Biofilm
