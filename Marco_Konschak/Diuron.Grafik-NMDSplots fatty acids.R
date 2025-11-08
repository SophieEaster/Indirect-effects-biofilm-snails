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

setwd("C:/Users/Marco Konschak/Desktop/Promotion/Methoden (Labor)/Fatty acids/GC in Landau/R-Skript für FA-Auswertung/Skripte für Studien/Grazer study")

Biofilm<-data.table(read.table("2019-10-27-Grazer-FAs of biofilm.txt"))
Biofilm[, Name := as.character(Name)]
Biofilm[, RUNINFO := as.character(RUNINFO)]
Biofilm[, Run := as.factor(Run)]
Biofilm[, Treatment := as.factor(Treatment)]
Biofilm[, Replicate := as.character(Replicate)]
Biofilm[, FA_mass := as.numeric(FA_mass)]
Biofilm[, µg_FA_per_mg := as.numeric(µg_FA_per_mg)]

Biofilm<-Biofilm[order(Treatment),]
Biofilm[,Name := gsub(" ", "", Name)]

Names<-unique(Biofilm$Name)
Names
RUNINFOs<-Biofilm$RUNINFO[Biofilm$Name=="16:0"]
Runs<-Biofilm$Run[Biofilm$Name=="16:0"]
Treatments<-Biofilm$Treatment[Biofilm$Name=="16:0"]
Replicates<-Biofilm$Replicate[Biofilm$Name=="16:0"]

`C15:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="15:0"]
`C16:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="16:0"]
`C16:1n-7` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="16:1n-7"]
`C17:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="17:0"]
`C18:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="18:0"]
`C20:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="20:0"]
`C18:1n-7` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="18:1n-7"]
`C18:1n-9c`<-Biofilm$µg_FA_per_mg[Biofilm$Name=="18:1n-9c"]
`C18:2n-6c`<-Biofilm$µg_FA_per_mg[Biofilm$Name=="18:2n-6c"]
`C18:3n-3` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="18:3n-3"]
`C18:3n-6` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="18:3n-6"]
`C20:1n-9` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="20:1n-9"]
`C20:2n-6` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="20:2n-6"]
`C20:3n-3` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="20:3n-3"]
`C20:4n-6` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="20:4n-6"]
`C20:5n-3` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="20:5n-3"]
`C22:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="22:0"]
`C22:6n-3` <-Biofilm$µg_FA_per_mg[Biofilm$Name=="22:6n-3"]
`C23:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="23:0"]
`C24:0`    <-Biofilm$µg_FA_per_mg[Biofilm$Name=="24:0"]

Biofilm2<-data.table(data.frame(RUNINFO=RUNINFOs,Run=Runs,Treatment=Treatments,Replicate=Replicates,
                                `C15:0`=`C15:0`,`C16:0`=`C16:0`,`C16:1n-7`=`C16:1n-7`,`C17:0`=`C17:0`,`C18:0`=`C18:0`,
                                `C18:1n-7`=`C18:1n-7`,`C18:1n-9c`=`C18:1n-9c`,`C18:2n-6c`=`C18:2n-6c`,`C18:3n-3`=`C18:3n-3`,
                                `C18:3n-6`=`C18:3n-6`,`C20:0`=`C20:0`,`C20:1n-9`=`C20:1n-9`,`C20:2n-6`=`C20:2n-6`,`C20:3n-3`=`C20:3n-3`,
                                `C20:4n-6`=`C20:4n-6`,`C20:5n-3`=`C20:5n-3`,`C22:0`=`C22:0`,`C22:6n-3`=`C22:6n-3`,`C23:0`=`C23:0`,`C24:0`=`C24:0`))


Biofilm2_ohneTR<-Biofilm2[,c(-1,-2,-3,-4)]


`X_C15:0`    <-Biofilm$FA_mass[Biofilm$Name=="15:0"]
`X_C16:0`    <-Biofilm$FA_mass[Biofilm$Name=="16:0"]
`X_C16:1n-7` <-Biofilm$FA_mass[Biofilm$Name=="16:1n-7"]
`X_C17:0`    <-Biofilm$FA_mass[Biofilm$Name=="17:0"]
`X_C18:0`    <-Biofilm$FA_mass[Biofilm$Name=="18:0"]
`X_C18:1n-7` <-Biofilm$FA_mass[Biofilm$Name=="18:1n-7"]
`X_C18:1n-9c`<-Biofilm$FA_mass[Biofilm$Name=="18:1n-9c"]
`X_C18:2n-6c`<-Biofilm$FA_mass[Biofilm$Name=="18:2n-6c"]
`X_C18:3n-3` <-Biofilm$FA_mass[Biofilm$Name=="18:3n-3"]
`X_C18:3n-6` <-Biofilm$FA_mass[Biofilm$Name=="18:3n-6"]
`X_C20:0`    <-Biofilm$FA_mass[Biofilm$Name=="20:0"]
`X_C20:1n-9` <-Biofilm$FA_mass[Biofilm$Name=="20:1n-9"]
`X_C20:2n-6` <-Biofilm$FA_mass[Biofilm$Name=="20:2n-6"]
`X_C20:3n-3` <-Biofilm$FA_mass[Biofilm$Name=="20:3n-3"]
`X_C20:4n-6` <-Biofilm$FA_mass[Biofilm$Name=="20:4n-6"]
`X_C20:5n-3` <-Biofilm$FA_mass[Biofilm$Name=="20:5n-3"]
`X_C22:0`    <-Biofilm$FA_mass[Biofilm$Name=="22:0"]
`X_C22:6n-3` <-Biofilm$FA_mass[Biofilm$Name=="22:6n-3"]
`X_C23:0`    <-Biofilm$FA_mass[Biofilm$Name=="23:0"]
`X_C24:0`    <-Biofilm$FA_mass[Biofilm$Name=="24:0"]


total_FAs<-data.table(data.frame(RUNINFO=RUNINFOs,Run=Runs,Treatment=Treatments,Replicate=Replicates,
                                 `C15:0`=`X_C15:0`,`C16:0`=`X_C16:0`,`C16:1n-7`=`X_C16:1n-7`,`C17:0`=`X_C17:0`,`C18:0`=`X_C18:0`,`C18:1n-7`=`X_C18:1n-7`,
                                 `C18:1n-9c`=`X_C18:1n-9c`,`C18:2n-6c`=`X_C18:2n-6c`, `C18:3n-3`=`X_C18:3n-3`,`C18:3n-6`=`X_C18:3n-6`,`C20:0`=`X_C20:0`,
                                 `C20:1n-9`=`X_C20:1n-9`,`C20:2n-6`=`X_C20:2n-6`,`C20:3n-3`=`X_C20:3n-3`,
                                 `C20:4n-6`=`X_C20:4n-6`,`C20:5n-3`=`X_C20:5n-3`,`C22:0`=`X_C22:0`,
                                 `C22:6n-3`=`X_C22:6n-3`,`C23:0`=`X_C23:0`,`C24:0`=`X_C24:0`))


###relativ
total_FAs_ohneTR<-total_FAs[,c(-1,-2,-3,-4)]
total_FAs_rel_ohneTR<-total_FAs_ohneTR/rowSums(total_FAs_ohneTR)*100
rowSums(total_FAs_rel_ohneTR)
total_FAs_rel<-cbind(total_FAs[,c(1,2,3,4)],total_FAs_rel_ohneTR)

Parameters<-data.table(read_excel("2019-03-07-Grazer-Environm. parameters.xlsx"))
Parameters[, Run := as.factor(Run)]
Parameters[, Treatment := as.factor(Treatment)]
Parameters[, Replicate := as.character(Replicate)]
Parameters[, ID := as.factor(ID)]
Parameters[, Concentration := as.numeric(Concentration)]
Parameters[, pH := as.numeric(pH)]
Parameters[, O2 := as.numeric(O2)]
Parameters[, Conductivity := as.numeric(Conductivity)]
Parameters[, Velocity := as.numeric(Velocity)]
Parameters[, Illuminance := as.numeric(Illuminance)]
Parameters[, PhotoActRad := as.numeric(PhotoActRad)]

#### Absolute 

Parameters<-Parameters[,c(-1,-2,-3)]
Biofilm2$ID<- with(Biofilm2, paste0(Run, Treatment, Replicate))

Biofilm2_par<-merge(Biofilm2, Parameters[, with = F],  by = "ID")
Biofilm2_par<-Biofilm2_par[order(Treatment),]

Parameters_mult<-subset(Biofilm2_par,select=c(Concentration,pH,O2,Conductivity,Velocity,Illuminance,PhotoActRad))
PhotoActRads<-Parameters_mult$PhotoActRad
Velocities<-Parameters_mult$Velocity
pHs<-Parameters_mult$pH
O2s<-Parameters_mult$O2
Conductivities<-Parameters_mult$Conductivity

Biofilm_MULT<-subset(Biofilm2_par,select=-c(ID,RUNINFO,Run, Treatment,Replicate,Concentration,pH,O2,Conductivity,Velocity,Illuminance,PhotoActRad))

FAs_mult<-sqrt(Biofilm_MULT)

#Selection of dimensions and Stress level (based and data transformation and resemblance metric)
FA_NMDS<-metaMDS(FAs_mult,k=2, autotransform = FALSE ,distance = "bray")
FA_NMDS$stress


#Fits environmental factors and vectors
colnames(FAs_mult)<-c("15:0","16:0","16:1n-7","17:0","18:0","18:1n-7","18:1n-9", "18:2n-6", "18:3n-3",  "18:3n-6",
                      "20:0","20:1n-9" ,"20:2n-6","20:3n-3","20:4n-6","20:5n-3","22:0","22:6n-3","23:0","24:0")
FA_variables.SAFAs<-envfit(FA_NMDS, subset(FAs_mult,select = c("15:0","16:0","17:0","18:0","20:0","22:0","23:0","24:0"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.MUFAs<-envfit(FA_NMDS, subset(FAs_mult,select = c("16:1n-7","18:1n-7","18:1n-9","20:1n-9"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.PUFAs<-envfit(FA_NMDS, subset(FAs_mult,select = c("18:2n-6","18:3n-3","18:3n-6","20:2n-6","20:3n-3","20:4n-6","20:5n-3","22:6n-3"))
                           ,permutations = 999, na.rm = TRUE)

# NMDS Plot

data.scores = as.data.frame(scores(FA_NMDS))
data.scores$Treatment = Treatments

SAFAs_coord_cont = as.data.frame(scores(FA_variables.SAFAs, "vectors")) * ordiArrowMul(FA_variables.SAFAs)
MUFAs_coord_cont = as.data.frame(scores(FA_variables.MUFAs, "vectors")) * ordiArrowMul(FA_variables.MUFAs)
PUFAs_coord_cont = as.data.frame(scores(FA_variables.PUFAs, "vectors")) * ordiArrowMul(FA_variables.PUFAs)

absol_PERI<-ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, size = 7, alpha = 0.5,shape=c(rep(21,2),rep(22,3),rep(24,3),rep(21,3),rep(22,3),rep(24,3)),
             fill=c(rep("white",8),rep("grey50",9)), color="black",stroke = 1.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = SAFAs_coord_cont, size =2, alpha = 0.5, colour = "darkorange") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = MUFAs_coord_cont, size =2, alpha = 0.5, colour = "deepskyblue3") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = PUFAs_coord_cont, size =2, alpha = 0.5, colour = "chartreuse3") +
  geom_text_repel(data = rbind(SAFAs_coord_cont,MUFAs_coord_cont,PUFAs_coord_cont), aes(x = NMDS1, y = NMDS2),
                  direction= "both", segment.size = 1, alpha = 1, size =5, colour = "black", fontface = "bold",
                  label = row.names(rbind(SAFAs_coord_cont,MUFAs_coord_cont,PUFAs_coord_cont)),
                  arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "first")) +
  theme(axis.title = element_text(size = 20, colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  annotate(geom="text", x=-0.6, y=0.2, label="a",color="black",size = 10)


#### Relative


total_FAs_rel$ID<- with(total_FAs_rel, paste0(Run, Treatment, Replicate))

total_FAs_rel_par<-merge(total_FAs_rel, Parameters[, with = F],  by = "ID")
total_FAs_rel_par<-total_FAs_rel_par[order(Treatment),]

total_FAs_rel_MULT<-subset(total_FAs_rel_par,select=-c(ID,RUNINFO,Run, Treatment,Replicate,Concentration,pH,O2,Conductivity,Velocity,Illuminance,PhotoActRad))
total_FAs_rel_MULT<-total_FAs_rel_MULT/100
relFAs_mult<-sqrt(total_FAs_rel_MULT)

# NMDS Plot

#Selection of dimensions and Stress level (based and data transformation and resemblance metric)
relFA_NMDS<-metaMDS(relFAs_mult,k=2, autotransform = FALSE ,distance = "bray")
relFA_NMDS$stress

#Fits environmental factors and vectors
colnames(relFAs_mult)<-c("15:0","16:0","16:1n-7","17:0","18:0","18:1n-7","18:1n-9", "18:2n-6", "18:3n-3",  "18:3n-6",
                         "20:0","20:1n-9" ,"20:2n-6","20:3n-3","20:4n-6","20:5n-3","22:0","22:6n-3","23:0","24:0")
relFA_variables.SAFAs<-envfit(relFA_NMDS, subset(relFAs_mult,select = c("15:0","16:0","17:0","18:0","20:0","22:0","23:0","24:0"))
                              ,permutations = 999, na.rm = TRUE)
relFA_variables.SAFAs

relFA_variables.MUFAs<-envfit(relFA_NMDS, subset(relFAs_mult,select = c("16:1n-7","18:1n-7","18:1n-9","20:1n-9"))
                              ,permutations = 999, na.rm = TRUE)
relFA_variables.MUFAs

relFA_variables.PUFAs<-envfit(relFA_NMDS, subset(relFAs_mult,select = c("18:2n-6","18:3n-3","18:3n-6","20:2n-6","20:3n-3","20:4n-6","20:5n-3","22:6n-3"))
                              ,permutations = 999, na.rm = TRUE)
relFA_variables.PUFAs

data.scores = as.data.frame(scores(relFA_NMDS))
data.scores$Treatment = Treatments

relSAFAs_coord_cont = as.data.frame(scores(relFA_variables.SAFAs, "vectors")) * ordiArrowMul(relFA_variables.SAFAs)
relMUFAs_coord_cont = as.data.frame(scores(relFA_variables.MUFAs, "vectors")) * ordiArrowMul(relFA_variables.MUFAs)
relPUFAs_coord_cont = as.data.frame(scores(relFA_variables.PUFAs, "vectors")) * ordiArrowMul(relFA_variables.PUFAs)

rel_PERI<-ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, size = 7, alpha = 0.5,shape=c(rep(21,2),rep(22,3),rep(24,3),rep(21,3),rep(22,3),rep(24,3)),
             fill=c(rep("white",8),rep("grey50",9)), color="black",stroke = 1.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = relSAFAs_coord_cont, size =2, alpha = 0.5, colour = "darkorange") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = relMUFAs_coord_cont, size =2, alpha = 0.5, colour = "deepskyblue3") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = relPUFAs_coord_cont, size =2, alpha = 0.5, colour = "chartreuse3") +
  geom_text_repel(data = rbind(relSAFAs_coord_cont,relMUFAs_coord_cont,relPUFAs_coord_cont), aes(x = NMDS1, y = NMDS2),
                  direction= "both", segment.size = 1,alpha = 1, size =5, colour = "black", fontface = "bold",
                  label = row.names(rbind(relSAFAs_coord_cont,relMUFAs_coord_cont,relPUFAs_coord_cont)),
                  arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "first")) +
  theme(axis.title = element_text(size = 20, colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  annotate(geom="text", x=-0.25, y=0.2, label="b",color="black",size = 10)



setwd("C:/Users/Marco Konschak/Desktop/Promotion/Methoden (Labor)/Fatty acids/GC in Landau/R-Skript für FA-Auswertung/Skripte für Studien/Grazer study")

Snail<-data.table(read.table("2019-10-27-Grazer-FAs of snails.txt"))
Snail[, Name := as.character(Name)]
Snail[, RUNINFO := as.character(RUNINFO)]
Snail[, Treatment := as.factor(Treatment)]
Snail[, Replicate := as.character(Replicate)]
Snail[, FA_mass := as.numeric(FA_mass)]
Snail[, µg_FA_per_mg := as.numeric(µg_FA_per_mg)]

Snail<-Snail[order(Treatment),]
Snail[,Name := gsub(" ", "", Name)]

Names<-unique(Snail$Name)
Names
RUNINFOs<-Snail$RUNINFO[Snail$Name=="16:0"]
Treatments<-Snail$Treatment[Snail$Name=="16:0"]
Replicates<-Snail$Replicate[Snail$Name=="16:0"]

`C12:0`    <-Snail$µg_FA_per_mg[Snail$Name=="12:0"]
`C14:0`    <-Snail$µg_FA_per_mg[Snail$Name=="14:0"]
`C15:0`    <-Snail$µg_FA_per_mg[Snail$Name=="15:0"]
`C16:0`    <-Snail$µg_FA_per_mg[Snail$Name=="16:0"]
`C16:1n-7` <-Snail$µg_FA_per_mg[Snail$Name=="16:1n-7"]
`C17:0`    <-Snail$µg_FA_per_mg[Snail$Name=="17:0"]
`C18:0`    <-Snail$µg_FA_per_mg[Snail$Name=="18:0"]
`C18:1n-7` <-Snail$µg_FA_per_mg[Snail$Name=="18:1n-7"]
`C18:1n-9c`<-Snail$µg_FA_per_mg[Snail$Name=="18:1n-9c"]
`C18:2n-6c`<-Snail$µg_FA_per_mg[Snail$Name=="18:2n-6c"]
`C18:3n-3` <-Snail$µg_FA_per_mg[Snail$Name=="18:3n-3"]
`C18:3n-6` <-Snail$µg_FA_per_mg[Snail$Name=="18:3n-6"]
`C20:0`    <-Snail$µg_FA_per_mg[Snail$Name=="20:0"]
`C20:1n-9` <-Snail$µg_FA_per_mg[Snail$Name=="20:1n-9"]
`C20:2n-6` <-Snail$µg_FA_per_mg[Snail$Name=="20:2n-6"]
`C20:3n-3` <-Snail$µg_FA_per_mg[Snail$Name=="20:3n-3"]
`C20:3n-6` <-Snail$µg_FA_per_mg[Snail$Name=="20:3n-6"]
`C20:4n-6` <-Snail$µg_FA_per_mg[Snail$Name=="20:4n-6"]
`C20:5n-3` <-Snail$µg_FA_per_mg[Snail$Name=="20:5n-3"]
`C22:0`    <-Snail$µg_FA_per_mg[Snail$Name=="22:0"]
`C22:1n-9` <-Snail$µg_FA_per_mg[Snail$Name=="22:1n-9"]
`C22:6n-3` <-Snail$µg_FA_per_mg[Snail$Name=="22:6n-3"]
`C23:0`    <-Snail$µg_FA_per_mg[Snail$Name=="23:0"]
`C24:0`    <-Snail$µg_FA_per_mg[Snail$Name=="24:0"]

Snail2<-data.table(data.frame(Treatment=Treatments,Replicate=Replicates,`C12:0`=`C12:0`,`C14:0`=`C14:0`,
                              `C15:0`=`C15:0`,`C16:0`=`C16:0`,`C16:1n-7`=`C16:1n-7`,`C17:0`=`C17:0`,
                              `C18:0`=`C18:0`,`C18:1n-7`=`C18:1n-7`,`C18:1n-9c`=`C18:1n-9c`,
                              `C18:2n-6c`=`C18:2n-6c`,`C18:3n-3`=`C18:3n-3`,`C18:3n-6`=`C18:3n-6`,
                              `C20:0`=`C20:0`,`C20:1n-9`=`C20:1n-9`,`C20:2n-6`=`C20:2n-6`,
                              `C20:3n-3`=`C20:3n-3`,`C20:3n-6`=`C20:3n-6`,`C20:4n-6`=`C20:4n-6`,
                              `C20:5n-3`=`C20:5n-3`,`C22:0`=`C22:0`,`C22:1n-9`=`C22:1n-9`,
                              `C22:6n-3`=`C22:6n-3`,`C23:0`=`C23:0`,`C24:0`=`C24:0`))

Snail2_ohneTR<-Snail2[,c(-1,-2)]


`X_C12:0`    <-Snail$FA_mass[Snail$Name=="12:0"]
`X_C14:0`    <-Snail$FA_mass[Snail$Name=="14:0"]
`X_C15:0`    <-Snail$FA_mass[Snail$Name=="15:0"]
`X_C16:0`    <-Snail$FA_mass[Snail$Name=="16:0"]
`X_C16:1n-7` <-Snail$FA_mass[Snail$Name=="16:1n-7"]
`X_C17:0`    <-Snail$FA_mass[Snail$Name=="17:0"]
`X_C18:0`    <-Snail$FA_mass[Snail$Name=="18:0"]
`X_C18:1n-7` <-Snail$FA_mass[Snail$Name=="18:1n-7"]
`X_C18:1n-9c`<-Snail$FA_mass[Snail$Name=="18:1n-9c"]
`X_C18:2n-6c`<-Snail$FA_mass[Snail$Name=="18:2n-6c"]
`X_C18:3n-3` <-Snail$FA_mass[Snail$Name=="18:3n-3"]
`X_C18:3n-6` <-Snail$FA_mass[Snail$Name=="18:3n-6"]
`X_C20:0`    <-Snail$FA_mass[Snail$Name=="20:0"]
`X_C20:1n-9` <-Snail$FA_mass[Snail$Name=="20:1n-9"]
`X_C20:2n-6` <-Snail$FA_mass[Snail$Name=="20:2n-6"]
`X_C20:3n-3` <-Snail$FA_mass[Snail$Name=="20:3n-3"]
`X_C20:3n-6` <-Snail$FA_mass[Snail$Name=="20:3n-6"]
`X_C20:4n-6` <-Snail$FA_mass[Snail$Name=="20:4n-6"]
`X_C20:5n-3` <-Snail$FA_mass[Snail$Name=="20:5n-3"]
`X_C22:0`    <-Snail$FA_mass[Snail$Name=="22:0"]
`X_C22:1n-9` <-Snail$FA_mass[Snail$Name=="22:1n-9"]
`X_C22:6n-3` <-Snail$FA_mass[Snail$Name=="22:6n-3"]
`X_C23:0`    <-Snail$FA_mass[Snail$Name=="23:0"]
`X_C24:0`    <-Snail$FA_mass[Snail$Name=="24:0"]


total_FAs<-data.table(data.frame(Treatment=Treatments,Replicate=Replicates,`C12:0`=`X_C12:0`,`C14:0`=`X_C14:0`,
                                 `C15:0`=`X_C15:0`,`C16:0`=`X_C16:0`,`C16:1n-7`=`X_C16:1n-7`,`C17:0`=`X_C17:0`,
                                 `C18:0`=`X_C18:0`,`C18:1n-7`=`X_C18:1n-7`,`C18:1n-9c`=`X_C18:1n-9c`,
                                 `C18:2n-6c`=`X_C18:2n-6c`,`C18:3n-3`=`X_C18:3n-3`,`C18:3n-6`=`X_C18:3n-6`,
                                 `C20:0`=`X_C20:0`,`C20:1n-9`=`X_C20:1n-9`,`C20:2n-6`=`X_C20:2n-6`,
                                 `C20:3n-3`=`X_C20:3n-3`,`C20:3n-6`=`X_C20:3n-6`,`C20:4n-6`=`X_C20:4n-6`,
                                 `C20:5n-3`=`X_C20:5n-3`,`C22:0`=`X_C22:0`,`C22:1n-9`=`X_C22:1n-9`,
                                 `C22:6n-3`=`X_C22:6n-3`,`C23:0`=`X_C23:0`,`C24:0`=`X_C24:0`))


###relativ
total_FAs_ohneTR<-total_FAs[,c(-1,-2)]
total_FAs_rel_ohneTR<-total_FAs_ohneTR/rowSums(total_FAs_ohneTR)*100
rowSums(total_FAs_rel_ohneTR)
total_FAs_rel<-cbind(total_FAs[,c(1,2)],total_FAs_rel_ohneTR)
Snail2_ohneTR<-Snail2[,c(-1,-2)]


#### Absolute

Snail2_mult<-sqrt(Snail2_ohneTR)

#Selection of dimensions and Stress level (based and data transformation and resemblance metric)
FA_NMDS<-metaMDS(Snail2_mult,k=2, autotransform = FALSE ,distance = "bray")
FA_NMDS$stress

#Fits environmental factors and vectors
colnames(Snail2_mult)<-c("12:0","14:0","15:0","16:0","16:1n-7","17:0","18:0","18:1n-7","18:1n-9", "18:2n-6", "18:3n-3",  "18:3n-6",
                         "20:0","20:1n-9" ,"20:2n-6","20:3n-3","20:3n-6","20:4n-6","20:5n-3","22:0","22:1n-9","22:6n-3","23:0","24:0")


FA_variables.SAFAs<-envfit(FA_NMDS, subset(Snail2_mult,select = c("12:0","14:0","15:0","16:0","17:0","18:0","20:0","22:0","23:0","24:0"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.MUFAs<-envfit(FA_NMDS, subset(Snail2_mult,select = c("16:1n-7","18:1n-7","18:1n-9","20:1n-9","22:1n-9"))
                           ,permutations = 999, na.rm = TRUE)

FA_variables.PUFAs<-envfit(FA_NMDS, subset(Snail2_mult,select = c("18:2n-6","18:3n-3","18:3n-6","20:2n-6","20:3n-3","20:3n-6","20:4n-6","20:5n-3","22:6n-3"))
                           ,permutations = 999, na.rm = TRUE)

# NMDS Plot

data.scores = as.data.frame(scores(FA_NMDS))
data.scores$Treatment = Treatments

SAFAs_coord_cont = as.data.frame(scores(FA_variables.SAFAs, "vectors")) * ordiArrowMul(FA_variables.SAFAs)
MUFAs_coord_cont = as.data.frame(scores(FA_variables.MUFAs, "vectors")) * ordiArrowMul(FA_variables.MUFAs)
PUFAs_coord_cont = as.data.frame(scores(FA_variables.PUFAs, "vectors")) * ordiArrowMul(FA_variables.PUFAs)

absol_SNAIL<-ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, size = 7, alpha = 0.5,shape=c(rep(21,6),rep(22,6)),
             fill=c(rep("white",6),rep("grey50",6)), color="black",stroke = 1.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = SAFAs_coord_cont, size =2, alpha = 0.5, colour = "darkorange") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = MUFAs_coord_cont, size =2, alpha = 0.5, colour = "deepskyblue3") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = PUFAs_coord_cont, size =2, alpha = 0.5, colour = "chartreuse3") +
  geom_text_repel(data = rbind(SAFAs_coord_cont,MUFAs_coord_cont,PUFAs_coord_cont), aes(x = NMDS1, y = NMDS2),
                  direction= "both", segment.size = 1,alpha = 1, size =5, colour = "black", fontface = "bold",
                  label = row.names(rbind(SAFAs_coord_cont,MUFAs_coord_cont,PUFAs_coord_cont)),
                  arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "first")) +
  theme(axis.title = element_text(size = 20, colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  annotate(geom="text", x=-0.25, y=0.2, label="c",color="black",size = 10)


#### Relative
relSnail_mult<-total_FAs_rel_ohneTR/100
relSnail_mult<-sqrt(relSnail_mult)

# NMDS Plot

#Selection of dimensions and Stress level (based and data transformation and resemblance metric)
relFA_NMDS<-metaMDS(relSnail_mult,k=2, autotransform = FALSE ,distance = "bray")
relFA_NMDS$stress

#Fits environmental factors and vectors
colnames(relSnail_mult)<-c("12:0","14:0","15:0","16:0","16:1n-7","17:0","18:0","18:1n-7","18:1n-9", "18:2n-6", "18:3n-3",  "18:3n-6",
                           "20:0","20:1n-9" ,"20:2n-6","20:3n-3","20:3n-6","20:4n-6","20:5n-3","22:0","22:1n-9","22:6n-3","23:0","24:0")


relFA_variables.SAFAs<-envfit(relFA_NMDS, subset(relSnail_mult,select = c("12:0","14:0","15:0","16:0","17:0","18:0","20:0","22:0","23:0","24:0"))
                              ,permutations = 999, na.rm = TRUE)

relFA_variables.MUFAs<-envfit(relFA_NMDS, subset(relSnail_mult,select = c("16:1n-7","18:1n-7","18:1n-9","20:1n-9","22:1n-9"))
                              ,permutations = 999, na.rm = TRUE)

relFA_variables.PUFAs<-envfit(relFA_NMDS, subset(relSnail_mult,select = c("18:2n-6","18:3n-3","18:3n-6","20:2n-6","20:3n-3","20:3n-6","20:4n-6","20:5n-3","22:6n-3"))
                              ,permutations = 999, na.rm = TRUE)


data.scores = as.data.frame(scores(relFA_NMDS))
data.scores$Treatment = Treatments

relSAFAs_coord_cont = as.data.frame(scores(relFA_variables.SAFAs, "vectors")) * ordiArrowMul(relFA_variables.SAFAs)
relMUFAs_coord_cont = as.data.frame(scores(relFA_variables.MUFAs, "vectors")) * ordiArrowMul(relFA_variables.MUFAs)
relPUFAs_coord_cont = as.data.frame(scores(relFA_variables.PUFAs, "vectors")) * ordiArrowMul(relFA_variables.PUFAs)

rel_SNAIL<-ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, size = 7, alpha = 0.5,shape=c(rep(21,6),rep(22,6)),
             fill=c(rep("white",6),rep("grey50",6)), color="black",stroke = 1.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = relSAFAs_coord_cont, size =2, alpha = 0.5, colour = "darkorange") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = relMUFAs_coord_cont, size =2, alpha = 0.5, colour = "deepskyblue3") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),data = relPUFAs_coord_cont, size =2, alpha = 0.5, colour = "chartreuse3") +
  geom_text_repel(data = rbind(relSAFAs_coord_cont,relMUFAs_coord_cont,relPUFAs_coord_cont), aes(x = NMDS1, y = NMDS2),
                  direction= "both", segment.size = 1,alpha = 1, size =5, colour = "black", fontface = "bold",
                  label = row.names(rbind(relSAFAs_coord_cont,relMUFAs_coord_cont,relPUFAs_coord_cont)),
                  arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "first")) +
  theme(axis.title = element_text(size = 20, colour = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "black"),
        axis.ticks = element_blank(), axis.text = element_blank()) +
  annotate(geom="text", x=-0.3, y=0.2, label="d",color="black",size = 10)


grid.arrange(absol_PERI, rel_PERI,absol_SNAIL,rel_SNAIL, ncol = 2)



