# NMDS for Snail + Biofilm (NL)FA
# 1. NLFA Snail absolute
# packages
require(vegan)
require(multcomp)
require(car)
require(data.table)
require(asbio)
require(plotrix)
library(readxl)
library(lmtest)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(dplyr)
library(tidyr)
library(shades)
library(ggpubr)


# insert data
# all NAs already omitted
NLFA_Snails_absolute <- read.table(file = 'C:\\Users\\sophi\\OneDrive\\Desktop\\Promotion\\Bioassays\\Schneckenversuch 2023\\FA_UFZ_Magdeburg\\Snail experiment\\NLFA_Snails_absolute.txt', 
                                   header = TRUE, sep = "\t")


View(NLFA_Snails_absolute)

#change to wide format
absolute_wide <- pivot_wider(NLFA_Snails_absolute, names_from=fa, values_from=fa_conc)

NMDS_absolute <- metaMDS(absolute_wide, k=2, autotransform = FALSE ,distance = "bray", na.rm = TRUE)
NMDS_absolute$stress








FA_variables.SAFAs<-envfit(FA_NMDS, subset(NLFA_Snails_absolute,select = c("C08:0", "C10:0", "C11:0", "C12:0", "C13:0", "C14:0", "C15:0","C16:0","C17:0","C18:0","C20:0", "C21:0", "C22:0","C24:0"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.MUFAs<-envfit(FA_NMDS, subset(NLFA_Snails_absolute,select = c("C14:1", "C15:1", "C16:1", "C17:1", "C18:1","C20:1","C21:1", "C22:1", "C24:1"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.PUFAs<-envfit(FA_NMDS, subset(NLFA_Snails_absolute,select = c("C18:2(n-6t)","C18:2(n-6c)", "C18:3(n-6)","C18:3(n-3)","C18:4(n-3)", "C20:2(n-6)","C20:3(n-6)","C20:4(n-6)","C20:5(n-3)","C22:5(n-3)", "C22:6(n-3)"))
                           ,permutations = 999, na.rm = TRUE)



data.scores = as.data.frame(scores(FA_NMDS)$sites)
data.scores$Treatment = Treatments


SAFAs_coord_cont = as.data.frame(scores(FA_variables.SAFAs, "vectors")) * ordiArrowMul(FA_variables.SAFAs)
MUFAs_coord_cont = as.data.frame(scores(FA_variables.MUFAs, "vectors")) * ordiArrowMul(FA_variables.MUFAs)
PUFAs_coord_cont = as.data.frame(scores(FA_variables.PUFAs, "vectors")) * ordiArrowMul(FA_variables.PUFAs)



#Selection of dimensions and Stress level (based and data transformation and resemblance metric)
library(vegan)
FA_NMDS <- metaMDS(data2_mult, k=2, autotransform = FALSE ,distance = "bray", na.rm = TRUE)
FA_NMDS$stress

# stress = 0.05443741


#Fits environmental factors and vectors
colnames(data2_mult)<-c("12:0", "13:0", "14:0", "14:1", "15:0", "15:1", "16:0","16:1n-7","17:0", "17:1", "18:0","18:1n-7","18:1n-9c", "18:2n-6", "18:3n-3",  "18:3n-6",
                        "20:0","20:1n-9" ,"20:2n-6","20:3n-3","20:3n-6","20:4n-6","20:5n-3", "21:0", "22:0","22:1n-9","22:6n-3","24:0", "24:1")


FA_variables.SAFAs<-envfit(FA_NMDS, subset(data2_mult,select = c("12:0", "13:0", "14:0", "15:0","16:0","17:0","18:0","20:0", "21:0", "22:0","24:0"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.MUFAs<-envfit(FA_NMDS, subset(data2_mult,select = c("14:1", "15:1", "16:1n-7", "17:1", "18:1n-7","18:1n-9c","20:1n-9","22:1n-9", "24:1"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.PUFAs<-envfit(FA_NMDS, subset(data2_mult,select = c("18:2n-6","18:3n-3","18:3n-6","20:2n-6","20:3n-3","20:3n-6","20:4n-6","20:5n-3","22:6n-3"))
                           ,permutations = 999, na.rm = TRUE)



data.scores = as.data.frame(scores(FA_NMDS)$sites)
data.scores$Treatment = Treatments


SAFAs_coord_cont = as.data.frame(scores(FA_variables.SAFAs, "vectors")) * ordiArrowMul(FA_variables.SAFAs)
MUFAs_coord_cont = as.data.frame(scores(FA_variables.MUFAs, "vectors")) * ordiArrowMul(FA_variables.MUFAs)
PUFAs_coord_cont = as.data.frame(scores(FA_variables.PUFAs, "vectors")) * ordiArrowMul(FA_variables.PUFAs)
#------------------------------------------------------------------------

plot(FA_NMDS)
plot(FA_variables.SAFAs)

absol_SNAIL = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, aes(colour = Treatment), size = 5, alpha = 0.5) +
  stat_conf_ellipse(data=data.scores, aes(colour = Treatment, fill=Treatment), alpha=0.3, geom="polygon") +
  
  
  #geom_point(data = data.scores, size = 7, alpha = 0.5,shape=c(rep(21,2),rep(22,3),rep(24,3),rep(21,3),rep(22,3),rep(24,3)),
  #fill=c(rep("white",5),rep("grey50",5),rep("grey10",5),rep("grey90",5)), color="black",stroke = 1.5) +
  # stat_conf_ellipse(geom = "polygon", aes(color = Treatment),alpha = 0.5) +
  #scale_colour_manual(values = c("orange", "steelblue", "chartreuse3", "violet")) +
  #geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
  # data = SAFAs_coord_cont, linewidth =1, alpha = 0.5, colour = "darkorange") +
  #geom_text(data = SAFAs_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", size=4, 
  #label = row.names(SAFAs_coord_cont)) +
  #geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
# data = MUFAs_coord_cont, linewidth =1, alpha = 0.5, colour = "deepskyblue3") +
#geom_text(data = MUFAs_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", size=4, 
#label = row.names(MUFAs_coord_cont)) +
#geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
#             data = PUFAs_coord_cont, linewidth =1, alpha = 0.5, colour = "black") +
#  geom_text(data = PUFAs_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", size=3, 
#            label = row.names(PUFAs_coord_cont)) +
theme_bw() +
  theme(axis.title = element_text(size = 11, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        #axis.ticks = element_blank(), 
        #axis.text = element_blank(), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 11, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 10, colour = "grey30"), 
        plot.margin = margin(t = 0.2, b = 0.2, l = 0.2, r = 0.2, unit = "cm"),
        panel.grid.major = element_line(colour = "grey95", size = 0.4),
        panel.grid.minor = element_line(colour = "grey95", size = 0.2)) +
  annotate("text", x = 0.4, y = -0.37,
           label=paste("stress:",FA_NMDS$stress %>% round(2)), size = 3)



absol_SNAIL 