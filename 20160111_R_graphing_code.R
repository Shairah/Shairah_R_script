

#generate random numbers to plot
set.seed(23)
x<-rnorm(100, mean=20, sd=3)
x#rnmorm is to give 100 of them?
y<-(x^2 + runif(100)) + rnorm(100, mean=20, sd=15)
y

group<-ifelse(x>20, "A", "B")

df <-data.frame(x,y, group)
df

head(df) #give the first 6 rows of the data

plot(x)  #give scatterplot
plot(x, y)
plot(y~x) #y is dependent var, x is independent var
#what's the difference betw histogram vs barplot
hist(x) #frequency of each level
barplot(x) #give the actual value

nd <-c(red=15, blue=10, green=22) #c mean concatenate
barplot(nd)

boxplot(x) #give median, with 25th and 75th quartile (for box) , whisker is the 2x of the data range of the quartile
#dot is outlier
#if dont want outlier, use boxplot(x, outline=F)
boxplot(x,y)
boxplot(df$x~df$group) ####need clarification

pairs(df)
z<- cbind(x,y) 
pairs(z)

#do subset
subA<-subset(df,group=="A")
plot(subA$x, subA$y)

############
#SCATTERPLOT#
plot(x,y,type="n")
#type="n" create plot without anything on it
#type ="l" create line graph
#type ="o" create line overlay (with point to the line) 
#type ="h" create histogram like

points(subA$x, subA$y, add=T) #hmmm

#add=T is used to add datapoint to the pre-existing graph

#labelling the graph
plot(x,y, main="X vs Y", xlab="x", ylab="y", pch=18, col="blue" )
plot(x,y, main="X vs Y", xlab="x", ylab="y", type="n" )
text(x,y,group,pos=4)

#pch is symbol for the point, each number represent different shape, diamond, dot, etc
#cex is size of the point label text (character expansion), default is 1.anything less than that give smaller text
#pos is position of the text label to the point (1,2,3,4 up, bottom, etc etc)
#lty is the line type

#Ctrl-C copy as bitmap create image taking the whole space
#Ctrl-W copy as metafile create clear/transparent background

#EXERCISE#
grow <-  read.csv( "~/SCHOOL/R_Graphing/growth.csv" )
head(grow)
x<-grow$age
y<-grow$length
plot(x,y)
plot(grow$age,grow$length)
plot(grow$length~grow$age)
boxplot(grow$length~grow$age)
before2008<-subset(grow, year>="2008")
head(before2008)
boxplot(before2008$length~before2008$age, main="length at age after yr 2007", cex.main="2", xlab="age", ylab="length")


##################################
#base package#

#quartz is to specify the width and height of the plot windows, can standardized size
plot(0)
quartz(w=4,h=4) #title???
#par is just giving the graphical parameter function
#mfrow and mfcol
#layout can also be used instead of par

par(mfrow=c(2,2))
plot(rnorm(100))
plot(rnorm(30))
plot(rnorm(70))
plot(rnorm(50))

m <-matrix(1:4, 2)
layout(m)
plot(rnorm(100))
plot(rnorm(30))
plot(rnorm(70))
plot(rnorm(50))

#layout.show will give the preview
ly<-layout(matrix(c(1,2,3,3), nrow=2, byrow=TRUE))
layout.show(ly)

#oma is the same like mar, but now we are talkng about the whole figure area
#c(bottom, left, top, right)
mtext("X AXIS", side=1, line=2, outer=TRUE) #device region is whole window?
      #label of title, side is position bottom, outer TRUE, put the label considering the the whole margin.

par(mar=c(0,0,0,0), oma=c(6,6,1,1), mfrow=c(2,3))
plot(rnorm(100))
plot(rnorm(30))
plot(rnorm(70))
plot(rnorm(50))
plot(rnorm(100))
plot(rnorm(30))

#xaxt is a axis type
par(mar=c(0,0,0,0), oma=c(6,6,1,1), mfrow=c(2,3))
plot(rnorm(50), xaxt='n') #xaxt='n' will not give the x axis label
plot(rnorm(50), xaxt='n', yaxt='n')
mtext(side=1, 'XAXIS', line=3)
mtext(side=1, 'XAXIS', line=3, outer=TRUE)
#axis this function can be used to add new axis?


##EXERCISE##
lob <-  read.csv( "~/SCHOOL/R_Graphing/lob.csv" )
head(lob)
subME<-subset(lob, STATE=="ME")
subMA<-subset(lob,STATE=="MA")
subCT<-subset(lob, STATE=="CT")
subNY<-subset(lob,STATE=="NY")

ly<-layout(matrix(c(1,2,3,4), nrow=2, byrow=TRUE))
layout.show(ly)
quartz(w=4,h=4)##questionable
par(mfrow=c(2,2), mar=c(0,0,0,0) , oma=c(6,6,2,2))
#par just setting the plot window
#par(mfrow) is the same with layout(matrix)

plot(subME$YEAR,subME$MT, log="y", xaxt='n', type="l")
plot(subMA$YEAR,subMA$MT, log="y", yaxt='n', xaxt='n', type="l")
plot(subCT$YEAR,subCT$MT, log="y")
plot(subNY$YEAR,subNY$MT, log="y", yaxt='n')

dev.off()
xaxt=FALSE #need to be inside the function argument




##################
#adding & modifying components of your graphs
#to close the active windows, use dev.off() 
plot(x,y,ylim=c(-200,1000), xlim=c(-5,30))
abline(h=0)
abline(v=0)

abline(reg)
reg=lm(y~x)


set.seed(23)
w<-rnorm(20,mean=19,sd=4)
z<-(-8*w+runif(20)+rnorm(20,mean=600,sd=5))
plot(x,y)
points(w,z,type='l') 

df2=data.frame(w,z)
df2<-df2[order(w),]
lines(df2$w,df2$z,type='l') #type='l' will just create a line, if you want to change the line type, have to specify lty="2"

plot(x,y)
lines(df2$w,df2$z,type='l', lty='3', col='purple', lwd=4, pch=18)
lines(df2$w,df2$z,type='l', col='blue', pch=15, cex=2)

plot(x,y)
curve(x^2+2*x+100,from=12,to=27,lty=4,lwd=3,add=TRUE)
somefunct <- function(x) 50*sin(x)+400
curve(somefunct,from=12,to=28,add=TRUE,n=500)
text(19,650,expression(km[2]))
text(15,600,expression(integral(x^2+y^2*partialdiff*x~~partialdiff*y, -infinity, 0), paste("This is random text"), adj=c(.25,.5)))


dev.off() #this function will reset the plot screen
# "\n command will give new line

plot(x,y)
text(x=23.5,625, labels=expression(y=x^2-3*x+100),srt=45)
curve(x^2-3*x+100, from =12, to =28, add=TRUE)

#at=c(15, 19,24) telling the plot the exact location on x axis where to put the text 
mtext(list("random", "word","placement"),side=1, at=c(15,19,24), line=4,cex=c(.5,.75,1.5)) #line=4 is the position at the bottom

##add legend##
plot(x,y,col='forestgreen', cex=1.75,lwd=2, pch=7)
abline(reg,lwd=5,lty=6,col='royalblue')
curve(somefunct,from=12,to=28,lty=4,lwd=3,add=TRUE,n=500,col='red')

legend('topleft',labels=expression("Data Points","Regression Line", y==50sin(x)+200),
col = c('forestgreen','royalblue','red'),pch=c(7,NA,NA),lty=c(0,6,4),lwd=c(NA,5,3),pt.cex=c(1.75,0,0),pt.lwd=c(2,0,0), horiz=TRUE))
#n is default for the number which is 100

####EXERCISE####
dev.off()
set.seed(23)
w<-rnorm(20,mean=19,sd=4)
z<-(-8*w+runif(20))+rnorm(20,mean=600,sd=5)

##############################
##COLORS

#colors()
#palette()

plot(x,y,pch=21,col="red")
plot(x,y,col=c(1:8))

##how to color by group of point
#factor in R has level
plot(x,y, col=as.factor(df$group))
plot(x,y, col=as.factor(df$group))
palette(rainbow(10))


plot(df$x,df$y)
plot(df$x,df$y, col=df$group)
plot(df$x,df$y, col=rank(df$x))

lisa.pal<-colorRampPalette(c("whitesmoke","hotpink", "green"))
palette(lisa.pal)

#ramppalette telling how many color that u want
#use R color brewer - use color developed by Cynthia Brewer
#RcolorBrewer

#combine RcolorBrewer + color ramp palette will give more color within that color brewer
library(RColorBrewer)
display.brewer.all()

#if you want to reset the palette to default use palette("default")

palette("default")

################################
##GGplot##

library('ggplot2')
library('plyr')

head(grow)
str(grow)
grow$age<-as.factor(grow$age)
grow$year<-as.factor(grow$year)

str(grow)

p_lw<-ggplot(data=grow,aes(x=length,y=weight))
p_lw+geom_blank()
p_lw+geom_point()
p_lw+geom_line()

a<-ggplot(grow,aes(length)) #need to specify the axis for histogram and density plot
a+geom_histogram()
a+geom_density()
#ddply kinda doing subsetting and summarise is doing some analyzing
mm<-ddply(grow,'age',summarise,mu=mean(length),sd=sd(length))
mm

barp<-ggplot(mm,aes(x=age,y=mu))+geom_bar(stat='identity')
#now want to add error bar
barp+geom_errorbar(aes(ymin=mu-2*sd, ymax=mu+2*sd),color='red')

##EXERCISE##
#make sure x axis is categorical data if we want to create plot in boxplot
p_lw2<-ggplot(data=grow,aes(x=age,y=length))
p_lw2+geom_boxplot()

ggplot(grow,aes(age,length))+stat_summary(fun.y=mean,geom="point")
ggplot(grow,aes(age,length))+stat_summary(fun.y=mean,geom="bar")

ggplot(grow,aes(age,length))+stat_summary(fun.y=mean,geom="point")+stat_summary(fun.data=mean_sdl,geom='errorbar',colour='red')

library('Himsc')
###EXERCISE###
ddply(grow,'age',summarise,med=median(length), high=quantile(length,0.7), low=quantile(length,0.3))

#got error

##use coord_flip() for transpose?
