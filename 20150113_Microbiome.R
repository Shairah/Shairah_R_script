library(ggplot2)
library(plyr)

SArare <- read.csv( "~/Desktop/20160107_midwest/SA-rarefaction1.csv" )

head(SArare)
str(SArare)

SArare2 <- ggplot(data=SArare,aes(x=length,y=weight))
SArare$numsampled <- as.factor(SArare$numsampled)
SArare$streamBS_t1 <- as.factor(SArare$streamBS_t1)
SArare$streamBS_t2  <- as.factor(SArare$streamBS_t2 ) 
SArare$streamBS_t3  <- as.factor(SArare$streamBS_t3 )

SArare$streamBSplus_t1  <- as.factor(SArare$streamBSplus_t1 )
SArare$streamBSplus_t2  <- as.factor(SArare$streamBSplus_t2 )
SArare$streamBSplus_t3  <- as.factor(SArare$streamBSplus_t3 )

SArare$gwUVBS_t1  <- as.factor(SArare$gwUVBS_t1 )
SArare$gwUVBS_t2  <- as.factor(SArare$gwUVBS_t2 )
SArare$gwUVBS_t3  <- as.factor(SArare$gwUVBS_t3 )

SArare$gwUVBSplus_t1  <- as.factor(SArare$gwUVBSplus_t1  )
SArare$gwUVBSplus_t2   <- as.factor(SArare$gwUVBSplus_t2 )
SArare$gwUVBSplus_t3   <- as.factor(SArare$gwUVBSplus_t3  )

p_rare<-ggplot(data=SArare,aes(x=numsampled,y=streamBS_t1, ))
p_rare+geom_line()

plot(SArare$numsampled, c(SArare$streamBS_t1, SArare$streamBSplus_t2, SArare$streamBSplus_t3 ),type="p")

####invsimpson###
Alpha <- read.csv( "~/Desktop/20160107_midwest/invsimp.csv" )
head(Alpha)
x<-time
y<-mean(subset(Alpha$invsimpson))
y
x<-rnorm(100, mean=20, sd=3)
x#rnmorm is to give 100 of them?
y<-(x^2 + runif(100)) + rnorm(100, mean=20, sd=15)
y

group<-ifelse(x>20, "A", "B")

df <-data.frame(x,y, group)
df


#############################################################################
#working code#
dev.off()
# Create Line Chart / rarefaction
Rare <- read.csv( "~/Desktop/Workbook11.csv", header=T )
head(Rare)
str(Rare)

# convert factor to numeric for convenience 
Rare$Treatment <- as.numeric(Rare$Treatment)
Rare$Treatment
ntrt <- max(Rare$Treatment)

# get the range for the x and y axis 
xrange <- range(Rare$num_sequence) 
yrange <- range(Rare$sobs) 

# set up the plot 
plot(xrange, yrange, type="n", xlab="no. of sequences",
     ylab="no. of OTU" )
colors <- rainbow(Rare$Treatment)
linetype <- c(1:ntrt) 
plotchar <- seq(18,18+ntrt,1)


# add lines 
for (i in 1:ntrt) { 
  Trt <- subset(Rare, Treatment==i) 
  lines(Trt$num_sequence, Trt$sobs,  type="l", lwd=5.0, lty=linetype[i],
       col=colors[i], pch=plotchar[i]) 
} 

# add a title and subtitle 
title("Rarefaction Curve by Treatment at different timepoint")
# add a legend 
legend(xrange[4], yrange[4], 1:ntrt, cex=5.0, col=colors, pch=plotchar, lty=linetype, 
         title="Treatment")
#issue with legend, and only 8 trt appear.. supposed to be 12

# Boxplot of indices by group 
ly<-layout(matrix(c(1,1,1,2,3,4), nrow=2, byrow=TRUE))
layout.show(ly)
colors()
#notimepoint
data7 <- read.csv( "~/Desktop/Workbook7.csv", header=T )
head(data7)
str(data7)
boxplot(invsimpson~group, data=data7, col=(c("brown", "yellow","gold", "darkgreen", "blue", "darkblue")), main="Diversity index", xaxt='n', ylab="Inv. Simpson, 1/D")
#timepoint1
data8 <- read.csv( "~/Desktop/Workbook8.csv", header=T )
head(data8)
boxplot(invsimpson~group, data=data8, col=(c("red4", "hotpink","lightblue","steelblue")), main="Diversity index", xlab="Group", ylab="Inv. Simpson, 1/D", ylim=c(0,23))

#timepoint2
data9 <- read.csv( "~/Desktop/Workbook9.csv", header=T )
head(data9)
boxplot(invsimpson~group, data=data9, col=(c("red4", "hotpink","lightblue","steelblue")), main="Diversity index", xlab="Group", ylab="Inv. Simpson, 1/D", ylim=c(0,23))
#timepoint3
data10 <- read.csv( "~/Desktop/Workbook10.csv", header=T )
head(data10)
boxplot(invsimpson~group, data=data10, col=(c("red4", "hotpink","lightblue","steelblue")), main="Diversity index", xlab="Group", ylab="Inv. Simpson, 1/D", ylim=c(0,23))
