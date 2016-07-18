#survival_analysis
library(km.ci) #Required package for KM curves
library(survival) #Package used for survival analysis functions
library(KMsurv)
library(help=KMsurv)
#create survival object data frame

####### Control 1 TankID: CT_1

id=paste("CT_1",c(1:50),sep="_")
Status=c(rep(0,29),rep(1,2),rep(0,19))
Time=c(rep(35,29),30,31,rep(35,19))

CT1_surv=data.frame(id,Time,Status)
CT1_surv

####### Control 2 TankID: CT_2

id=paste("CT_2",c(1:50),sep="_")
Status=rep(0,50)
Time=rep(35,50)

CT2_surv=data.frame(id,Time,Status)


####### Control 3 TankID: CT_3

id=paste("CT_3",c(1:50),sep="_")
Status=rep(0,50)
Time=rep(35,50)

CT3_surv=data.frame(id,Time,Status)

####### Treatment 1 TankID: TR_1

id=paste("TR_1",c(1:50),sep="_")
Status=c(rep(0,30),rep(1,20))
Time=c(rep(35,30),32,33,33,33,rep(34,7),rep(35,9))

TR1_surv=data.frame(id,Time,Status)

####### Treatment 2 TankID: TR_2

id=paste("TR_2",c(1:50),sep="_")
Status=c(rep(0,29),rep(1,21))
Time=c(rep(35,29),32,32,33,33,33,rep(34,7),rep(35,9))

TR2_surv=data.frame(id,Time,Status)

####### Treatment 3 TankID: TR_3

id=paste("TR_3",c(1:50),sep="_")
Status=c(rep(0,20),rep(1,30))
Time=c(rep(35,20),2,2,30,31,rep(32,6),rep(33,9),rep(34,10),35)

TR3_surv=data.frame(id,Time,Status)

################## All Data
Group=c(rep("CT",150),rep("TR",150))

FullSurvival=rbind(CT1_surv,CT2_surv,CT3_surv,TR1_surv,TR2_surv,TR3_surv)
FullSurvival$Group=Group

##survival analysis
TrFit <- survfit(Surv(Time, Status) ~ Group, data = FullSurvival)
plot(TrFit, col=c("coral1", "darkturquoise"), lty = 1:1, lwd=2, ylab="Percent of survival", xlab="Days of exogenous feeding", conf.int =TRUE) 
legend("bottomright", c("CT","TR"), col=c("coral1", "darkturquoise"), lty = 1:1, lwd =2)
#To add confidence intervals, add conf.int=TRUE after ylim=c(0,1)

summary(TrFit) #Summary data by fruit

#Is there a difference between Group variables?
#This is the log-rank test
TrFit2<-survdiff(Surv(Time, Status)~Group, data=FullSurvival)
TrFit2
#p < .01 implies yes.

     

####TRlength####
getwd
setwd("/Users/shairahar/Dropbox/RESEARCH/progress/DietTransition_Chapter_3")
TRmm <-read.csv("Transition_Length.csv", header=T )
head(TRmm)
dim(TRmm)
View(TRmm)

##Try to create bar plot showing Length by Group at each time using ggplot 

# Boxplots of Length by group at 3 sampling point
library(ggplot2)
=qplot(data = TRmm, 
      x = Group, fill=Group,  main="Length by Group at 3 sampling time",
      facets = ~Time)





#have some issue onhow to put label for X and Y axis, and Title
Dboxplot=ggplot(data = TRmm) + 
  geom_boxplot(aes(y=Length, x=Group, fill=Group, ylab="Length(mm)")) + 
  facet_wrap(~Time)
Dboxplot+labs(y="Length(mm)")
#some tips
#http://docs.ggplot2.org/0.9.2.1/labs.html


###One way anova treat time as cont var (Transition vs Control- not distinguish between treatment at each time)
LengthCTR<-lm(Length~Group+nTime-1, data=TRmm)
LengthCTR
summary(LengthCTR)
Lengthnotime<-lm(Length~Group, data=TRmm)

trmm14=TRmm[TRmm$nTime==14,]
trmm21=TRmm[TRmm$nTime==21,]
trmm35=TRmm[TRmm$nTime==35,]

t.test(Length~Group,trmm14)
t.test(Length~Group,trmm21)
t.test(Length~Group,trmm35)

summary(lm(Length~Group,trmm14))
summary(LengthCTR)
##working code###

#################

##trial and/or possibly failed code##

#do subset
subT14<-subset(TRmm,Time=="14dpf")
head(subT14)
tail(subT14)
subT21<-subset(TRmm,Time=="21dpf")
subT35<-subset(TRmm,Time=="35dpf")

dev.off()
ly<-layout(matrix(c(1,2,3), nrow=1, byrow=TRUE))
par(mfrow=c(1,3), mai = c(1, 1, 0.1, 0.1))
layout.show(ly)


layout.show(op)

#boxplot(responsevar~bywhat, data=datafile)

boxplot(Length~Group, data=subT14, col=(c("yellow", "brown")), cex.axis = "1.5", ylab="Length, mm", cex.lab="2", ylim=c(25,60))
boxplot(Length~Group, data=subT21, col=(c("yellow", "brown")), cex.axis="1.5", yaxt='n', ylim=c(25,60))
boxplot(Length~Group, data=subT35, col=(c("yellow","brown")), cex.axis="1.5", yaxt='n', ylim=c(25,60))
title(main="Length for Control and Transitioned group at T14, T21, T35", cex.main="3"), xlab="Groups", cex.lab="2"

##some code taken from website
###Example of codes
mydata <- data.frame(myGroup = c('a', 'b'), myX = c(1,1))
qplot(data = mydata, 
      x = myX, 
      facets = ~myGroup)
ggplot(data = mydata) + 
  geom_bar(aes(myX)) + 
  facet_wrap(~myGroup)

# observations (points) are overlayed and jittered
qplot(Group, Length, data=subT14, geom="boxplot", 
      fill=Group, main="Length by Group",
      xlab="", ylab="Length(mm)")



library(reshape)
dm <- melt(d, id.vars = "var8")
dm$var8 <- as.factor(dm$var8)

ggplot(dm, aes(x = var8, y = value, fill = var8)) + 
  theme_bw() +
  scale_fill_manual(values = colors, guide = FALSE) +
  geom_boxplot()+ facet_wrap(~ variable)


TRfish<-matrix(c(rep(9,300*5)),300,5)
View(TRfish)
colnames<-c("fishID","status","day","group","replicate")
colnames(TRfish)<-c("fishID","status","day","group","replicate")
head(TRfish)

TRfish[,1]<-seq(1:300)
head(TRfish)

TRfish[,4]<-c(rep("CT",150),rep("TR",150))
View(TRfish)
TRfish[,2]<-c(rep(1,974),rep(0,1276),rep(1,743),rep(0,1179))
View(data)

TRfish[,5]<-c(rep(1,50),rep(2,100),rep(3,150),rep(1,200),rep(2,250),rep(3,300))
#have some problem on how to denote the replicate


data[,2]<-c(rep(1,974),rep(0,1276),rep(1,743),rep(0,1179))
View(data)
data[,3]<-c(rep(1,643),rep(2,231),rep(3,35),rep(4,32),rep(5,19),rep(6,2),rep(7,1),rep(8,11),rep(9,0),rep(9,1276),    rep(1,224),rep(2,181),rep(3,170),rep(4,80),rep(5,52),rep(6,26),rep(7,5),rep(8,4),rep(9,1),rep(9,1179)
            
