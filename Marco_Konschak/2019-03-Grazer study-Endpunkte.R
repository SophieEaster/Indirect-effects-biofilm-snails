library(vegan)
library(multcomp)
library(car)
library(data.table)
library(asbio)
library(plotrix)
library(graphics)
library(drc)
library(lsr)
library(readxl)
library(nlme)

daten_21<-data.table(read_excel("2019-09-Grazer study-Frate over 21 days.xlsx"))
daten_21[, Day := as.numeric(Day)]
daten_21[, Treatment := as.factor(Treatment)]
daten_21[, Replicate := as.character(Replicate)]
daten_21[, Subject := as.numeric(Subject)]
daten_21[, Frate := as.numeric(Frate)]
daten_21[, Area := as.numeric(Area)]

daten_21<-subset(daten_21,! Frate %in% NA)
daten_21<-within(daten_21,log_Frate<-log10(daten_21$Frate/mean(daten_21$Frate)+0.1))
daten_21<-within(daten_21,log_Area<-log10(daten_21$Area/mean(daten_21$Area)+0.1))

#===================
#### Frate in mg/day
#===================

boxplot(daten_21$Frate[daten_21$Day=="7"]~daten_21$Treatment[daten_21$Day=="7"],range=1.5,main="Day 7")
boxplot(daten_21$Frate[daten_21$Day=="14"]~daten_21$Treatment[daten_21$Day=="14"],range=1.5,main="Day 14")
boxplot(daten_21$Frate[daten_21$Day=="21"]~daten_21$Treatment[daten_21$Day=="21"],range=1.5,main="Day 21")

daten_21<-daten_21[order(Subject),]
hist(daten_21$Frate) # right skewed
hist(daten_21$log_Frate)

#AV_daten_21<-aov(log_Frate~factor(Day)*Treatment+Error(factor(Subject)),data=daten_21)

intercept      <-gls(log_Frate ~ 1, data = daten_21, method = "ML",na.action = na.omit)
randomintercept<-lme(log_Frate ~ 1, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)
timeFE         <-lme(log_Frate ~ 1 + Day, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)
timeRS         <-lme(log_Frate ~ 1 + Day, random = ~Day|Subject, data = daten_21, method = "ML",na.action = na.omit)
anova(intercept,randomintercept,timeFE,timeRS)

Frate_baseline     <-lme(log_Frate ~ 1, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

Frate_Model_t      <-lme(log_Frate ~ 1 + Day, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

`Frate_Model_t+Tox`<-lme(log_Frate ~ Day + Treatment, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

`Frate_Model_t*Tox` <-lme(log_Frate ~ Day * Treatment, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

anova(Frate_baseline,Frate_Model_t,`Frate_Model_t+Tox`,`Frate_Model_t*Tox`)
Model<-`Frate_Model_t*Tox`
summary(Model)
intervals(Model, 0.95)

Levene.Model<-lm(resid(Model) ~ daten_21$Subject)
anova(Levene.Model)
plot(Model)
plot(fitted(Model),resid(Model))

qqnorm(Model, ~ranef(., level=1))
qqPlot(Model$residuals)
hist(residuals(Model,type = "normalized"))
shapiro.test(residuals(Model,type = "normalized"))

#### Pairwise comparisons
daten_7<-subset(daten_21, Day %in% "7")
tapply(daten_7$Frate,daten_7$Treatment,shapiro.test)

par(mfrow=c(1,2))
qqPlot(daten_7$Frate[daten_7$Treatment=="C"])
qqPlot(daten_7$Frate[daten_7$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_7$Frate,daten_7$Treatment)

wilcox.test(daten_7$Frate[daten_7$Treatment=="T"],daten_7$Frate[daten_7$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)


daten_14<-subset(daten_21, Day %in% "14")
tapply(daten_14$Frate,daten_14$Treatment,shapiro.test)

par(mfrow=c(1,2))
qqPlot(daten_14$Frate[daten_14$Treatment=="C"])
qqPlot(daten_14$Frate[daten_14$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_14$Frate,daten_14$Treatment)

wilcox.test(daten_14$Frate[daten_14$Treatment=="T"],daten_14$Frate[daten_14$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)

daten_end<-subset(daten_21, Day %in% "21")
tapply(daten_end$Frate,daten_end$Treatment,shapiro.test)

par(mfrow=c(1,2))
qqPlot(daten_end$Frate[daten_end$Treatment=="C"])
qqPlot(daten_end$Frate[daten_end$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_end$Frate,daten_end$Treatment)

wilcox.test(daten_end$Frate[daten_end$Treatment=="T"],daten_end$Frate[daten_end$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)
# sig.

par(mfrow=c(1,2))
plot(daten_21$Day[daten_21$Treatment=="C"],daten_21$log_Frate[daten_21$Treatment=="C"],xlab="Day",ylab="log(Frate)",main="Control")
plot(daten_21$Day[daten_21$Treatment=="T"],daten_21$log_Frate[daten_21$Treatment=="T"],xlab="Day",ylab="log(Frate)",main="Treatment")
par(mfrow=c(1,1))
Control<-lm(daten_21$log_Frate[daten_21$Treatment=="C"]~daten_21$Day[daten_21$Treatment=="C"])
#plot(Control)
Treatment<-lm(daten_21$log_Frate[daten_21$Treatment=="T"]~daten_21$Day[daten_21$Treatment=="T"])
#plot(Treatment)

#### Grafik

medians.C<-tapply(daten_21$Frate[daten_21$Treatment=="C"],daten_21$Day[daten_21$Treatment=="C"],median)
medians.T<-tapply(daten_21$Frate[daten_21$Treatment=="T"],daten_21$Day[daten_21$Treatment=="T"],median)

lower.C<-tapply(daten_21$Frate[daten_21$Treatment=="C"],daten_21$Day[daten_21$Treatment=="C"],function(x) ci.median(x)$ci[2])
lower.T<-tapply(daten_21$Frate[daten_21$Treatment=="T"],daten_21$Day[daten_21$Treatment=="T"],function(x) ci.median(x)$ci[2])

upper.C<-tapply(daten_21$Frate[daten_21$Treatment=="C"],daten_21$Day[daten_21$Treatment=="C"],function(x) ci.median(x)$ci[3])
upper.T<-tapply(daten_21$Frate[daten_21$Treatment=="T"],daten_21$Day[daten_21$Treatment=="T"],function(x) ci.median(x)$ci[3])

# Control
x1<-c(1,2,3)
par(mar=c(4.5, 5.5, 2, 4) + 0.1,xpd=TRUE)
plotCI(x1,medians.C,ui=upper.C,li=lower.C,ylim=c(0,3),xlim=c(1,3.2),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=21,
       col="black",ylab=expression(""),xlab="Day",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
#lines(x1,medians.C,type="l",lty=1,col="black",lwd=1.5)
axis(1,at=c(1.1125,2.1125,3.1125),labels=c("7","14","21"),tick=TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=c(0.0,0.5,1.0,1.5,2.0,2.5,3.0),labels=c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
lines(c(0.9,3.4),c(0,0),type="l",lty=1,lwd=1.5)
mtext("Feeding rate in mg/day",at=1.5,side=2,line=4,cex=1.5)

# Water
x2<-c(1.25,2.25,3.25)
plotCI(x2,medians.T,ui=upper.T,li=lower.T,ylim=c(0,3),las=1,yaxs="i",pch=22,col="black",pt.bg="grey50",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5,add=TRUE)
#lines(x2,medians.T,type="l",lty=2,col="grey50",lwd=1.5)
lines(c(2.95,3.3),c(2.65,2.65),type="l",lty=1,lwd=1.5)
text(3.125,2.75,"*",cex=2)


# Grafik für Verteidigung

# Control
x1<-c(1,2,3)
par(mar=c(4.5, 5.5, 2, 4) + 0.1,xpd=TRUE)
plotCI(x1,medians.C,ui=upper.C,li=lower.C,ylim=c(0,3),xlim=c(1,3.2),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",
       ylab=expression(""),xlab="Tag",pch=21,col="black", pt.bg ="green",cex=2,lwd=2,cex.lab=1.5)

axis(1,at=c(1.1125,2.1125,3.1125),labels=c("7","14","21"),tick=TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=c(0.0,0.5,1.0,1.5,2.0,2.5,3.0),labels=c("0.0","0.5","1.0","1.5","2.0","2.5","3.0"),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
lines(c(0.9,3.4),c(0,0),type="l",lty=1,lwd=1.5)
mtext("Periphyton-Konsum in mg/d",at=1.5,side=2,line=4,cex=1.5)

# Water
x2<-c(1.25,2.25,3.25)
plotCI(x2,medians.T,ui=upper.T,li=lower.T,ylim=c(0,3),las=1,yaxs="i",pch=21,col="black", pt.bg ="orange",cex=2,lwd=2,add=TRUE)

# 6 x 9


#========
#### Area
#========

par(mar=c(4.5, 5.5, 2, 4) + 0.1,xpd=FALSE)
linear_reg.model<-lm(Frate~Area,data=daten_21)
#plot(linear_reg.model)
plot(daten_21$Area,daten_21$Frate,xlim=c(0,22),ylim=c(0,12),xaxt = "n",yaxt = "n",bty ="n",las=1,yaxs="i",pch=21,
     xlab=expression("Grazed area in cm"^"2"),ylab="Feeding rate in mg/d",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
axis(1,at=c(0,2,4,6,8,10,12,14,16,18,20,22),labels=c(0,2,4,6,8,10,12,14,16,18,20,22),tick=TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=c(0,2,4,6,8,10,12),labels=c(0,2,4,6,8,10,12),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
abline(linear_reg.model,lwd=1.5)
newx<-seq(min(daten_21$Area),25,length.out=1000)
conf_bands<-predict(linear_reg.model, newdata=data.frame("Area"=newx), interval="confidence",level = 0.95)
lines(newx,conf_bands[,2], col="blue", lty=2,lwd=1.5)
lines(newx, conf_bands[,3], col="blue", lty=2,lwd=1.5)
lines(c(-1,23),c(0,0),type="l",lty=1,lwd=1.5)
# 8 x 6
summary(lm(daten_21$Frate~daten_21$Area))

boxplot(daten_21$Area[daten_21$Day=="7"]~daten_21$Treatment[daten_21$Day=="7"],range=1.5,main="Day 7")
boxplot(daten_21$Area[daten_21$Day=="14"]~daten_21$Treatment[daten_21$Day=="14"],range=1.5,main="Day 14")
boxplot(daten_21$Area[daten_21$Day=="21"]~daten_21$Treatment[daten_21$Day=="21"],range=1.5,main="Day 21")

daten_21<-daten_21[order(Subject),]
hist(daten_21$Area)
hist(daten_21$log_Area)


intercept      <-gls(log_Area ~ 1, data = daten_21, method = "ML",na.action = na.omit)
randomintercept<-lme(log_Area ~ 1, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)
timeFE         <-lme(log_Area ~ 1 + Day, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)
timeRS         <-lme(log_Area ~ 1 + Day, random = ~Day|Subject, data = daten_21, method = "ML",na.action = na.omit)
anova(intercept,randomintercept,timeFE,timeRS)

Area_baseline     <-lme(log_Area ~ 1, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

Area_Model_t      <-lme(log_Area ~ 1 + Day, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

`Area_Model_t+Tox`<-lme(log_Area ~ Day + Treatment, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

`Area_Model_t*Tox` <-lme(log_Area ~ Day * Treatment, random = ~1|Subject, data = daten_21, method = "ML",na.action = na.omit)

anova(Area_baseline,Area_Model_t,`Area_Model_t+Tox`,`Area_Model_t*Tox`)
Model<-`Area_Model_t*Tox`
summary(Model)
intervals(Model, 0.95)

Levene.Model<-lm(resid(Model) ~ daten_21$Subject)
anova(Levene.Model)
plot(Model)
plot(fitted(Model),resid(Model))

qqnorm(Model, ~ranef(., level=1))
qqPlot(Model$residuals)
hist(residuals(Model,type = "normalized"))
shapiro.test(residuals(Model,type = "normalized"))


#### Pairwise comparisons 
daten_7<-subset(daten_21, Day %in% "7")
tapply(daten_7$Area,daten_7$Treatment,shapiro.test)

par(mfrow=c(1,2))
qqPlot(daten_7$Area[daten_7$Treatment=="C"])
qqPlot(daten_7$Area[daten_7$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_7$Area,daten_7$Treatment)

wilcox.test(daten_7$Area[daten_7$Treatment=="T"],daten_7$Area[daten_7$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)


daten_14<-subset(daten_21, Day %in% "14")
tapply(daten_14$Area,daten_14$Treatment,shapiro.test)

par(mfrow=c(1,2))
qqPlot(daten_14$Area[daten_14$Treatment=="C"])
qqPlot(daten_14$Area[daten_14$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_14$Area,daten_14$Treatment)

wilcox.test(daten_14$Area[daten_14$Treatment=="T"],daten_14$Area[daten_14$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)

daten_end<-subset(daten_21, Day %in% "21")
tapply(daten_end$Area,daten_end$Treatment,shapiro.test)

par(mfrow=c(1,2))
qqPlot(daten_end$Area[daten_end$Treatment=="C"])
qqPlot(daten_end$Area[daten_end$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_end$Area,daten_end$Treatment)

wilcox.test(daten_end$Area[daten_end$Treatment=="T"],daten_end$Area[daten_end$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)
# sig.

par(mfrow=c(1,2))
plot(daten_21$Day[daten_21$Treatment=="C"],daten_21$log_Area[daten_21$Treatment=="C"],xlab="Day",ylab="log(Area)",main="Control")
plot(daten_21$Day[daten_21$Treatment=="T"],daten_21$log_Area[daten_21$Treatment=="T"],xlab="Day",ylab="log(Area)",main="Treatment")
par(mfrow=c(1,1))
Control<-lm(daten_21$log_Area[daten_21$Treatment=="C"]~daten_21$Day[daten_21$Treatment=="C"])
#plot(Control)
Treatment<-lm(daten_21$log_Area[daten_21$Treatment=="T"]~daten_21$Day[daten_21$Treatment=="T"])
#plot(Treatment)


#### Grafik

medians.C<-tapply(daten_21$Area[daten_21$Treatment=="C"],daten_21$Day[daten_21$Treatment=="C"],median)
medians.T<-tapply(daten_21$Area[daten_21$Treatment=="T"],daten_21$Day[daten_21$Treatment=="T"],median)

lower.C<-tapply(daten_21$Area[daten_21$Treatment=="C"],daten_21$Day[daten_21$Treatment=="C"],function(x) ci.median(x)$ci[2])
lower.T<-tapply(daten_21$Area[daten_21$Treatment=="T"],daten_21$Day[daten_21$Treatment=="T"],function(x) ci.median(x)$ci[2])

upper.C<-tapply(daten_21$Area[daten_21$Treatment=="C"],daten_21$Day[daten_21$Treatment=="C"],function(x) ci.median(x)$ci[3])
upper.T<-tapply(daten_21$Area[daten_21$Treatment=="T"],daten_21$Day[daten_21$Treatment=="T"],function(x) ci.median(x)$ci[3])

# Control
x1<-c(1,2,3)
par(mar=c(4.5, 5.5, 2, 4) + 0.1,xpd=TRUE)
plotCI(x1,medians.C,ui=upper.C,li=lower.C,ylim=c(0,8),xlim=c(1,3.2),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=21,
       col="black",ylab=expression(""),xlab="Day",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
#lines(x1,medians.C,type="l",lty=1,col="black",lwd=1.5)
axis(1,at=c(1.1125,2.1125,3.1125),labels=c("7","14","21"),tick=TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=c(0,1,2,3,4,5,6,7,8),labels=c(0,1,2,3,4,5,6,7,8),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
lines(c(0.9,3.4),c(0,0),type="l",lty=1,lwd=1.5)
mtext(expression("Grazed area in cm"^"2"),at=4,side=2,line=4,cex=1.5)

# Water
x2<-c(1.25,2.25,3.25)
plotCI(x2,medians.T,ui=upper.T,li=lower.T,ylim=c(0,10),las=1,yaxs="i",pch=22,col="black",pt.bg="grey50",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5,add=TRUE)
#lines(x2,medians.T,type="l",lty=2,col="grey50",lwd=1.5)


#===============================================================================#

daten<-data.table(read_excel("2019-09-Grazer study-Endpunkte Bioassay.xlsx"))
daten[, Snail := as.character(Snail)]
daten[, Treatment := as.factor(Treatment)]
daten[, `growth/d_dry` := as.numeric(`growth/d_dry`)]
daten[, `Mean_Frate` := as.numeric(`Mean_Frate`)]
daten[, `Mean_Area` := as.numeric(`Mean_Area`)]
daten[, `Consumed_diuron` := as.numeric(`Consumed_diuron`)]

##############----dry weight----##############

daten_d.weight<-subset(daten,select=c(Treatment,Snail,`growth/d_dry`))
daten_d.weight<-subset(daten_d.weight,! `growth/d_dry` %in% NA)

boxplot(`growth/d_dry`~Treatment,daten_d.weight,range=1.5,main="growth over 21 d")

tapply(daten_d.weight$`growth/d_dry`,daten_d.weight$Treatment,shapiro.test) # C sig.

par(mar=c(4.5, 5.5, 2, 4) + 0.1,xpd=FALSE)
par(mfrow=c(1,2))
qqPlot(daten_d.weight$`growth/d_dry`[daten_d.weight$Treatment=="C"])
qqPlot(daten_d.weight$`growth/d_dry`[daten_d.weight$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_d.weight$`growth/d_dry`,daten_d.weight$Treatment)

wilcox.test(daten_d.weight$`growth/d_dry`[daten_d.weight$Treatment=="T"],daten_d.weight$`growth/d_dry`[daten_d.weight$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)

medians_d.weight<-tapply(daten_d.weight$`growth/d_dry`,daten_d.weight$Treatment,median)*1000
lower_d.weight<-tapply(daten_d.weight$`growth/d_dry`,daten_d.weight$Treatment,function(x) ci.median(x)$ci[2])*1000
upper_d.weight<-tapply(daten_d.weight$`growth/d_dry`,daten_d.weight$Treatment,function(x) ci.median(x)$ci[3])*1000

par(mar=c(5, 6.5, 2, 4) + 0.1,xpd=TRUE)
x<-c(1,2)
plotCI(x,medians_d.weight,ui=upper_d.weight,li=lower_d.weight,ylim=c(0,140),xlim=c(0,3),bty="n",las=1,yaxs="i",xaxt="n",yaxt = "n",ylab="",xlab="Diuron in µg/L",pch=16,cex.axis=1.5,cex.lab=1.5,lwd=1.5)
axis(side=1, at=c(1,2), labels=c("0","8"),las=1,cex.axis=1.5,cex.lab=1.5,lwd=1.5)
axis(2,at=c(0,20,40,60,80,100,120,140),labels=c(0,20,40,60,80,100,120,140),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
mtext("Dry weight in µg/day",at=70,side=2,line=4,cex=1.5)
lines(c(-0.15,3),c(0,0),type="l",lty=1,lwd=1.5)

##############----Frate----##############

daten_Frate<-subset(daten,select=c(Treatment,Snail,`Mean_Frate`))
daten_Frate<-subset(daten_Frate,! `Mean_Frate` %in% NA)

boxplot(`Mean_Frate`~Treatment,daten_Frate,range=1.5,main="Consumption per day")

tapply(daten_Frate$`Mean_Frate`,daten_Frate$Treatment,shapiro.test) # sig

par(mar=c(4.5, 5.5, 2, 4) + 0.1,xpd=FALSE)
par(mfrow=c(1,2))
qqPlot(daten_Frate$`Mean_Frate`[daten_Frate$Treatment=="C"])
qqPlot(daten_Frate$`Mean_Frate`[daten_Frate$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_Frate$`Mean_Frate`,daten_Frate$Treatment)

wilcox.test(daten_Frate$`Mean_Frate`[daten_Frate$Treatment=="T"],daten_Frate$`Mean_Frate`[daten_Frate$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)

medians_Frate<-tapply(daten_Frate$`Mean_Frate`,daten_Frate$Treatment,median)
lower_Frate<-tapply(daten_Frate$`Mean_Frate`,daten_Frate$Treatment,function(x) ci.median(x)$ci[2])
upper_Frate<-tapply(daten_Frate$`Mean_Frate`,daten_Frate$Treatment,function(x) ci.median(x)$ci[3])

par(mar=c(5, 6.5, 2, 4) + 0.1,xpd=TRUE)
x<-c(1,2)
plotCI(x,medians_Frate,ui=upper_Frate,li=lower_Frate,ylim=c(0,2),xlim=c(0,3),bty="n",las=1,yaxs="i",xaxt="n",yaxt = "n",ylab="",xlab="Diuron in µg/L",pch=16,cex.axis=1.5,cex.lab=1.5,lwd=1.5)
axis(side=1, at=c(1,2), labels=c("0","8"),las=1,cex.axis=1.5,cex.lab=1.5,lwd=1.5)
axis(2,at=c(0,0.5,1.0,1.5,2),labels=c("0.0","0.5","1.0","1.5","2.0"),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
mtext("Consumption in mg/day",at=1,side=2,line=4,cex=1.5)
lines(c(-0.15,3),c(0,0),type="l",lty=1,lwd=1.5)


##############----Grazed Area----##############

daten_Area<-subset(daten,select=c(Treatment,Snail,`Mean_Area`))
daten_Area<-subset(daten_Area,! `Mean_Area` %in% NA)

boxplot(`Mean_Area`~Treatment,daten_Area,range=1.5,main="Grazed Area per week")
par(mar=c(5, 6.5, 2, 4) + 0.1,xpd=F)
plot(daten_Area$Mean_Area,daten_Frate$`Mean_Frate`,xlab="Area",ylab="Frate")
abline(lm(daten_Frate$`Mean_Frate`~daten_Area$Mean_Area))
summary(lm(daten_Frate$`Mean_Frate`~daten_Area$Mean_Area))

tapply(daten_Area$`Mean_Area`,daten_Area$Treatment,shapiro.test) # sig.

par(mfrow=c(1,2))
qqPlot(daten_Area$`Mean_Area`[daten_Area$Treatment=="C"])
qqPlot(daten_Area$`Mean_Area`[daten_Area$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(daten_Area$`Mean_Area`,daten_Area$Treatment)

wilcox.test(daten_Area$`Mean_Area`[daten_Area$Treatment=="T"],daten_Area$`Mean_Area`[daten_Area$Treatment=="C"],alternative="two.sided",paired=FALSE,conf.level=1-0.05,conf.int=TRUE)

medians_Area<-tapply(daten_Area$`Mean_Area`,daten_Area$Treatment,median)
lower_Area<-tapply(daten_Area$`Mean_Area`,daten_Area$Treatment,function(x) ci.median(x)$ci[2])
upper_Area<-tapply(daten_Area$`Mean_Area`,daten_Area$Treatment,function(x) ci.median(x)$ci[3])

par(mar=c(5, 6.5, 2, 4) + 0.1,xpd=TRUE)
x<-c(1,2)
plotCI(x,medians_Area,ui=upper_Area,li=lower_Area,ylim=c(0,6),xlim=c(0,3),bty="n",las=1,yaxs="i",xaxt="n",yaxt = "n",ylab="",xlab="Diuron in µg/L",pch=16,cex.axis=1.5,cex.lab=1.5,lwd=1.5)
axis(side=1, at=c(1,2), labels=c("0","8"),las=1,cex.axis=1.5,cex.lab=1.5,lwd=1.5)
axis(2,at=c(0,1,2,3,4,5,6),labels=c(0,1,2,3,4,5,6),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
mtext("Grazed area in m^2 per week",at=3,side=2,line=4,cex=1.5)
lines(c(-0.15,3),c(0,0),type="l",lty=1,lwd=1.5)


##############----Diuron intake (µg during the entire study)----##############

daten_DIU<-subset(daten,select=c(Treatment,Snail,`Consumed_diuron`))
daten_DIU<-subset(daten_DIU,! `Consumed_diuron` %in% NA)

median(daten_DIU$Consumed_diuron)
ci.median(daten_DIU$Consumed_diuron)

##############----Mibis----##############


Mibi_daten<-data.table(read_excel("2019-03-Grazer study-Endpunkte Mibis.xlsx"))
Mibi_daten[, Run := as.factor(Run)]
Mibi_daten[, Treatment := as.factor(Treatment)]
Mibi_daten[, Replicate := as.character(Replicate)]
Mibi_daten[, ID := as.character(ID)]
Mibi_daten[, Photo_Light := as.numeric(Photo_Light)]
Mibi_daten[, Carotenoids := as.numeric(Carotenoids)]
Mibi_daten[, Chlorophyll_a := as.numeric(Chlorophyll_a)]
Mibi_daten[, Chlorophyll_b := as.numeric(Chlorophyll_b)]
Mibi_daten[, Chlorophyll_c := as.numeric(Chlorophyll_c)]
Mibi_daten[, Fluorescence := as.numeric(Fluorescence)]
Mibi_daten[, AFDM := as.numeric(AFDM)]

boxplot(Photo_Light~Treatment,Mibi_daten,range=1.5,main="Light")
boxplot(Carotenoids~Treatment,Mibi_daten,range=1.5,main="Carotenoids")
boxplot(Chlorophyll_a~Treatment,Mibi_daten,range=1.5,main="Chlorophyll_a")
boxplot(Chlorophyll_b~Treatment,Mibi_daten,range=1.5,main="Chlorophyll_b")
boxplot(Chlorophyll_c~Treatment,Mibi_daten,range=1.5,main="Chlorophyll_c")
boxplot(Fluorescence~Treatment,Mibi_daten,range=1.5,main="Fluorescence")
boxplot(AFDM~Treatment,Mibi_daten,range=1.5,main="AFDM")

par(xpd=FALSE)
plot(Mibi_daten$Fluorescence[Mibi_daten$Treatment=="C"],Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="C"],xlab="Fluorescence",ylab="Chl a",main="Control")
abline(lm(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="C"]~Mibi_daten$Fluorescence[Mibi_daten$Treatment=="C"]))
summary(lm(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="C"]~Mibi_daten$Fluorescence[Mibi_daten$Treatment=="C"]))

plot(Mibi_daten$Fluorescence[Mibi_daten$Treatment=="T"],Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="T"],xlab="Fluorescence",ylab="Chl a",main="Treatment")
abline(lm(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="T"]~Mibi_daten$Fluorescence[Mibi_daten$Treatment=="T"]))
summary(lm(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="T"]~Mibi_daten$Fluorescence[Mibi_daten$Treatment=="T"]))

AV_Photo_Light<-aov(Photo_Light~Treatment,Mibi_daten)
summary(AV_Photo_Light)

# Carotenoids
leveneTest(Mibi_daten$Carotenoids,Mibi_daten$Treatment)

AV_Carotenoids_Slopes<-aov(Carotenoids~Photo_Light*Treatment,Mibi_daten)
summary(AV_Carotenoids_Slopes) 
AV_Carotenoids<-aov(Carotenoids~Photo_Light+Treatment,Mibi_daten)
shapiro.test(AV_Carotenoids$residuals)
summary(AV_Carotenoids)
#plot(AV_Carotenoids)


# Chlorophyll_a
# t-Test / Wilcox
tapply(Mibi_daten$Chlorophyll_a,Mibi_daten$Treatment,shapiro.test)
par(mfrow=c(1,2))
qqPlot(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="C"])
qqPlot(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(Mibi_daten$Chlorophyll_a,Mibi_daten$Treatment)

t.test(Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="T"],Mibi_daten$Chlorophyll_a[Mibi_daten$Treatment=="C"],alternative="two.sided",paired=FALSE,var.equal=TRUE,conf.level=0.95)

#ANCOVA
leveneTest(Mibi_daten$Chlorophyll_a,Mibi_daten$Treatment)

AV_Chlorophyll_a_Slopes<-aov(Chlorophyll_a~Photo_Light*Treatment,Mibi_daten)
summary(AV_Chlorophyll_a_Slopes) 
AV_Chlorophyll_a<-aov(Chlorophyll_a~Photo_Light+Treatment,Mibi_daten)
summary(AV_Chlorophyll_a)
shapiro.test(AV_Chlorophyll_a$residuals)
#plot(AV_Chlorophyll_a)

# Two-way ANCOVA
leveneTest(Mibi_daten$Chlorophyll_a,interaction(Mibi_daten$Treatment,Mibi_daten$Run))
AV_Chlorophyll_a_Slopes<-aov(Chlorophyll_a~Photo_Light*Treatment*Run,Mibi_daten)
summary(AV_Chlorophyll_a_Slopes)
AV_Chlorophyll_a<-aov(Chlorophyll_a~Photo_Light+Treatment*Run,Mibi_daten)
shapiro.test(AV_Chlorophyll_a$residuals)
summary(AV_Chlorophyll_a)

# Chlorophyll_b
# t-Test / Wilcox
tapply(Mibi_daten$Chlorophyll_b,Mibi_daten$Treatment,shapiro.test)
par(mfrow=c(1,2))
qqPlot(Mibi_daten$Chlorophyll_b[Mibi_daten$Treatment=="C"])
qqPlot(Mibi_daten$Chlorophyll_b[Mibi_daten$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(Mibi_daten$Chlorophyll_b,Mibi_daten$Treatment)

t.test(Mibi_daten$Chlorophyll_b[Mibi_daten$Treatment=="T"],Mibi_daten$Chlorophyll_b[Mibi_daten$Treatment=="C"],alternative="two.sided",paired=FALSE,var.equal=TRUE,conf.level=0.95)

#ANCOVA
leveneTest(Mibi_daten$Chlorophyll_b,Mibi_daten$Treatment)

AV_Chlorophyll_b_Slopes<-aov(Chlorophyll_b~Photo_Light*Treatment,Mibi_daten)
summary(AV_Chlorophyll_b_Slopes) 
AV_Chlorophyll_b<-aov(Chlorophyll_b~Photo_Light+Treatment,Mibi_daten)
summary(AV_Chlorophyll_b)
shapiro.test(AV_Chlorophyll_b$residuals)
#plot(AV_Chlorophyll_b)

# Two-way ANCOVA
leveneTest(Mibi_daten$Chlorophyll_b,interaction(Mibi_daten$Treatment,Mibi_daten$Run))
AV_Chlorophyll_b_Slopes<-aov(Chlorophyll_b~Photo_Light*Treatment*Run,Mibi_daten)
summary(AV_Chlorophyll_b_Slopes)
AV_Chlorophyll_b<-aov(Chlorophyll_b~Photo_Light+Treatment*Run,Mibi_daten)
shapiro.test(AV_Chlorophyll_b$residuals)
summary(AV_Chlorophyll_b)

# Chlorophyll_c
# t-Test / Wilcox
tapply(Mibi_daten$Chlorophyll_c,Mibi_daten$Treatment,shapiro.test)
par(mfrow=c(1,2))
qqPlot(Mibi_daten$Chlorophyll_c[Mibi_daten$Treatment=="C"])
qqPlot(Mibi_daten$Chlorophyll_c[Mibi_daten$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(Mibi_daten$Chlorophyll_c,Mibi_daten$Treatment)

t.test(Mibi_daten$Chlorophyll_c[Mibi_daten$Treatment=="T"],Mibi_daten$Chlorophyll_c[Mibi_daten$Treatment=="C"],alternative="two.sided",paired=FALSE,var.equal=TRUE,conf.level=0.95)


#ANCOVA
leveneTest(Mibi_daten$Chlorophyll_c,Mibi_daten$Treatment)

AV_Chlorophyll_c_Slopes<-aov(Chlorophyll_c~Photo_Light*Treatment,Mibi_daten)
summary(AV_Chlorophyll_c_Slopes) 
AV_Chlorophyll_c<-aov(Chlorophyll_c~Photo_Light+Treatment,Mibi_daten)
summary(AV_Chlorophyll_c)
shapiro.test(AV_Chlorophyll_c$residuals)
#plot(AV_Chlorophyll_c)

# Two-way ANCOVA
leveneTest(Mibi_daten$Chlorophyll_c,interaction(Mibi_daten$Treatment,Mibi_daten$Run))
AV_Chlorophyll_c_Slopes<-aov(Chlorophyll_c~Photo_Light*Treatment*Run,Mibi_daten)
summary(AV_Chlorophyll_c_Slopes)
AV_Chlorophyll_c<-aov(Chlorophyll_c~Photo_Light+Treatment*Run,Mibi_daten)
summary(AV_Chlorophyll_c)
shapiro.test(AV_Chlorophyll_c$residuals)

# Fluorescence
leveneTest(Mibi_daten$Fluorescence,Mibi_daten$Treatment)

AV_Fluorescence_Slopes<-aov(Fluorescence~Photo_Light*Treatment,Mibi_daten)
summary(AV_Fluorescence_Slopes) 
AV_Fluorescence<-aov(Fluorescence~Photo_Light+Treatment,Mibi_daten)
summary(AV_Fluorescence)
shapiro.test(AV_Fluorescence$residuals)
#plot(AV_Fluorescence)

# AFDM
# t-Test / Wilcox
tapply(Mibi_daten$AFDM,Mibi_daten$Treatment,shapiro.test)
par(mfrow=c(1,2))
qqPlot(Mibi_daten$AFDM[Mibi_daten$Treatment=="C"])
qqPlot(Mibi_daten$AFDM[Mibi_daten$Treatment=="T"])
par(mfrow=c(1,1))

leveneTest(Mibi_daten$AFDM,Mibi_daten$Treatment)

t.test(Mibi_daten$AFDM[Mibi_daten$Treatment=="T"],Mibi_daten$AFDM[Mibi_daten$Treatment=="C"],alternative="two.sided",paired=FALSE,var.equal=TRUE,conf.level=0.95)

#ANCVOA
leveneTest(Mibi_daten$AFDM,Mibi_daten$Treatment)

AV_AFDM_Slopes<-aov(AFDM~Photo_Light*Treatment,Mibi_daten)
summary(AV_AFDM_Slopes) 
AV_AFDM<-aov(AFDM~Photo_Light+Treatment,Mibi_daten)
summary(AV_AFDM)
shapiro.test(AV_AFDM$residuals)
#plot(AV_AFDM)

# Two-way ANCOVA
leveneTest(Mibi_daten$AFDM,interaction(Mibi_daten$Treatment,Mibi_daten$Run))

AV_AFDM_Slopes<-aov(AFDM~Photo_Light*Treatment*Run,Mibi_daten)
summary(AV_AFDM_Slopes)
AV_AFDM<-aov(AFDM~Photo_Light+Treatment*Run,Mibi_daten)
summary(AV_AFDM)
shapiro.test(AV_AFDM$residuals)

medians_a<-tapply(Mibi_daten$Chlorophyll_a,Mibi_daten$Treatment,median)
lower_a<-tapply(Mibi_daten$Chlorophyll_a,Mibi_daten$Treatment,function(x) ci.median(x)$ci[2])
upper_a<-tapply(Mibi_daten$Chlorophyll_a,Mibi_daten$Treatment,function(x) ci.median(x)$ci[3])

medians_b<-tapply(Mibi_daten$Chlorophyll_b,Mibi_daten$Treatment,median)
lower_b<-tapply(Mibi_daten$Chlorophyll_b,Mibi_daten$Treatment,function(x) ci.median(x)$ci[2])
upper_b<-tapply(Mibi_daten$Chlorophyll_b,Mibi_daten$Treatment,function(x) ci.median(x)$ci[3])

medians_c<-tapply(Mibi_daten$Chlorophyll_c,Mibi_daten$Treatment,median)
lower_c<-tapply(Mibi_daten$Chlorophyll_c,Mibi_daten$Treatment,function(x) ci.median(x)$ci[2])
upper_c<-tapply(Mibi_daten$Chlorophyll_c,Mibi_daten$Treatment,function(x) ci.median(x)$ci[3])

medians_Carotenoids<-tapply(Mibi_daten$Carotenoids,Mibi_daten$Treatment,median)
lower_Carotenoids<-tapply(Mibi_daten$Carotenoids,Mibi_daten$Treatment,function(x) ci.median(x)$ci[2])
upper_Carotenoids<-tapply(Mibi_daten$Carotenoids,Mibi_daten$Treatment,function(x) ci.median(x)$ci[3])

medians_AFDM<-tapply(Mibi_daten$AFDM,Mibi_daten$Treatment,median)
lower_AFDM<-tapply(Mibi_daten$AFDM,Mibi_daten$Treatment,function(x) ci.median(x)$ci[2])
upper_AFDM<-tapply(Mibi_daten$AFDM,Mibi_daten$Treatment,function(x) ci.median(x)$ci[3])

medians_Fluo<-tapply(Mibi_daten$Fluorescence,Mibi_daten$Treatment,median)/1000000000
lower_Fluo<-tapply(Mibi_daten$Fluorescence,Mibi_daten$Treatment,function(x) ci.median(x)$ci[2])/1000000000
upper_Fluo<-tapply(Mibi_daten$Fluorescence,Mibi_daten$Treatment,function(x) ci.median(x)$ci[3])/1000000000

par(mfrow=c(2,1))
# AFDM
par(mar=c(2.5, 5.5, 2, 5.5) + 0.1,xpd=TRUE)
x4<-c(1,3)
plotCI(x4,medians_AFDM,ui=upper_AFDM,li=lower_AFDM,ylim=c(0,30),xlim=c(0,4),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=21,
       col="black",pt.bg=c("white","grey50"),ylab=expression(""),xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
axis(1,at=c(1,3),labels=c("",""),tick=TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=c(0,5,10,15,20,25,30),labels=c(0,5,10,15,20,25,30),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
mtext(expression(paste("AFDM in g/m"^"2")),at=15,side=2,line=3.5,cex=1.5)
lines(c(-0.2,4),c(0,0),type="l",lty=1,lwd=1.5)
text(0.15,30*0.95,"a",cex=2)

# Pigments
par(mar=c(4.5, 5.5, 0, 5.5) + 0.1,xpd=TRUE)
x1<-c(1,4)
plotCI(x1,medians_a,ui=upper_a,li=lower_a,ylim=c(0,350),xlim=c(0,6),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=22,
       col="black",pt.bg=c("white","grey50"),ylab=expression(""),xlab="Diuron concentration in µg/L",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
axis(1,at=c(1.5,4.5),labels=c("0","8"),tick=TRUE,cex.axis=1.5,lwd=1.5)
axis(2,at=c(0,50,100,150,200,250,300,350),labels=c(0,50,100,150,200,250,300,350),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
mtext(expression(paste("Chl ", italic("a")," in mg/m"^"2")),at=175,side=2,line=3.5,cex=1.5)
text(0.15,350*0.95,"b",cex=2)
par(new=TRUE)

x2<-c(1.5,4.5)
plotCI(x2,medians_b,ui=upper_b,li=lower_b,ylim=c(0,150),xlim=c(0,6),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=23,
       col="black",pt.bg=c("white","grey50"),ylab=expression(""),xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
axis(4,at=c(0,30,60,90,120,150),labels=c(0,30,60,90,120,150),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
par(new=TRUE)

x3<-c(2,5)
plotCI(x3,medians_c,ui=upper_c,li=lower_c,ylim=c(0,150),xlim=c(0,6),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=24,
       col="black",pt.bg=c("white","grey50"),ylab=expression(""),xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)


mtext(expression(paste("Chl ", italic("b"), " & ", italic("c")," in mg/m"^"2")),at=75,side=4,line=4,cex=1.5)
lines(c(-0.35,6.35),c(0,0),type="l",lty=1,lwd=1.5)

par(mfrow=c(1,1))

# 6 x 8

# x5<-c(2,5)
# plotCI(x5,medians_Fluo,ui=upper_Fluo,li=lower_Fluo,ylim=c(0,12),xlim=c(0,6),bty ="n",xaxt = "n",yaxt = "n",las=1,yaxs="i",pch=22,
#        col="black",pt.bg=c("white","grey50"),ylab=expression(""),xlab="",cex=1.5,cex.lab=1.5,cex.axis=1.5,lwd=1.5)
# axis(4,at=c(0,2,4,6,8,10,12),labels=c(0,expression(paste("2"^"9")),expression(paste("4"^"9")),
#                                       expression(paste("6"^"9")),expression(paste("8"^"9")),expression(paste("10"^"9")),
#                                       expression(paste("12"^"9"))),tick=TRUE,las=1,cex.axis=1.5,lwd=1.5)
# mtext(expression(paste("Fluorescence in units/m"^"2")),at=6,side=4,line=3.5,cex=1.5)
# lines(c(-0.35,6.35),c(0,0),type="l",lty=1,lwd=1.5)


