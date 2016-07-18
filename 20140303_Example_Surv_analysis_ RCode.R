setwd("//cstatnas.cstat.msu.edu/redirects/frisbymi/Desktop") #Set working directory
#Install packages by going to packages -> install packages -> USA(MI) -> select package
library(km.ci) #Required package for KM curves
library(survival) #Package used for survival analysis functions

#Read the data into a table in R
data <- read.csv("SurvivalFruitTemplate.csv", header=T)
data #display data to confirm its form

#Creation of survival object
fit <- survfit(Surv(DegreeDays, Status)~Fruit, data) 

#Plotting the survival curves by group
plot(fit, col=c("red", "blue", "green", "black"), lty = 1:1, 
ylab="Percent Without First Emergence", xlab="Degree Days", 
xlim=c(100,400), ylim=c(0,1))
#To add confidence intervals, add conf.int=TRUE after ylim=c(0,1)

#Legend
legend("topright", c("Cherry", "Honeycrisp", "GD-S", "Undetermined"), 
lty=1, col=c("red", "blue", "green", "black"))

#Title
title(main="Survival Curves Sorted by Fruit")

#Replotting with confidence intervals
Cherry <- survfit(Surv(DegreeDays,Status)~Fruit, conf.type = "plain", data)
plot(Cherry)

summary(fit) #Summary data by fruit

#Is there a difference between the straitification variables? (Fruits)
#This is the log-rank test
survdiff(Surv(DegreeDays, Status)~Fruit, data)
#p < .01 implies yes.

#Is there a difference between the straitification variables? (Fruits)
#This is the peto-modification of the Gehan-Wilcoxon test
survdiff(Surv(DegreeDays, Status)~Fruit, data, rho=1)
#p < .01 implies yes.