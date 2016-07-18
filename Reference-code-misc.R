## load vegan
require("vegan")

## load the Dune data
data(dune, dune.env)
dune
dune.env

## PCA of the Dune data
mod <- rda(dune, scale = TRUE)

## plot the PCA
plot(mod, scaling=3)

with(dune.env, levels(Use))
scl <- 3 ## scaling = 3

colvec <- c("red2", "green4", "mediumblue")

plot(mod, type = "n", scaling = scl)
with(dune.env, points(mod, display = "sites", col = colvec[Use],
                      scaling = scl, pch = 21, bg = colvec[Use]))

head(with(dune.env, colvec[Use]))
text(mod, display = "species", scaling = scl, cex = 0.8, col = "darkcyan")
with(dune.env, legend("topright", legend = levels(Use), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec))

##load ggplot
library(ggplot2)
##trial with ggplot
data(Orange)
Orange
qplot(age, circumference, data = Orange, geom = "line",
      colour = Tree,
      main = "How does orange tree circumference vary with age?")

#mtcars
data(mtcars)
mtcars
aggregate(mtcars$hp, by=list(carb=mtcars$carb, am=mtcars$am), mean)
aggregate(hp ~ carb + am, mtcars, mean) #alt formula
#last column is the number of carburettors (“carb”), which could be an good grouping variable. 
#Another is  whether the car is an automatic or not (“am”). 
#As our dependant variable we’ll use the horsepower (“hp”)

View(mtcars)
myData <- aggregate(mtcars$mpg,
                    by = list(cyl = mtcars$cyl, gears = mtcars$gear),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n = length(x)))
myData <- do.call(data.frame, myData)
myData$se <- myData$x.sd / sqrt(myData$x.n)

colnames(myData) <- c("cyl", "gears", "mean", "sd", "n", "se")

myData$names <- c(paste(myData$cyl, "cyl /",
                       myData$gears, " gear"))
#plot
par(mar = c(5, 6, 4, 5) + 0.1)

plotTop <- max(myData$mean) +
  myData[myData$mean == max(myData$mean), 6] * 3

barCenters <- barplot(height = myData$mean,
                      names.arg = myData$names,
                      beside = true, las = 2,
                      ylim = c(0, plotTop),
                      cex.names = 0.75, xaxt = "n",
                      main = "Mileage by No. Cylinders and No. Gears",
                      ylab = "Miles per Gallon",
                      border = "black", axes = TRUE)

# Specify the groupings. We use srt = 45 for a
# 45 degree string rotation
text(x = barCenters, y = par("usr")[3] - 1, srt = 45,
     adj = 1, labels = myData$names, xpd = TRUE)

segments(barCenters, myData$mean - myData$se * 2, barCenters,
         myData$mean + myData$se * 2, lwd = 1.5)

arrows(barCenters, myData$mean - myData$se * 2, barCenters,
       myData$mean + myData$se * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)
#grouped barplot
head(myData)
tapply(myData$mean, list(myData$cyl, myData$gears),
       function(x) c(x = x))
tabbedMeans <- tapply(myData$mean, list(myData$cyl,
                                        myData$gears),
                      function(x) c(x = x))
tabbedSE <- tapply(myData$se, list(myData$cyl,
                                   myData$gears),
                   function(x) c(x = x))

barCenters <- barplot(height = tabbedMeans,
                      beside = TRUE, las = 1,
                      ylim = c(0, plotTop),
                      cex.names = 0.75,
                      main = "Mileage by No. Cylinders and No. Gears",
                      ylab = "Miles per Gallon",
                      xlab = "No. Gears",
                      border = "black", axes = TRUE,
                      legend.text = TRUE,
                      args.legend = list(title = "No. Cylinders", 
                                         x = "topright",
                                         cex = .7))

segments(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
         tabbedMeans + tabbedSE * 2, lwd = 1.5)

arrows(barCenters, tabbedMeans - tabbedSE * 2, barCenters,
       tabbedMeans + tabbedSE * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

#trial
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

y <- rnorm(500, mean=1)
y <- matrix(y,100,5)
y.means <- apply(y,2,mean)
y.sd <- apply(y,2,sd)
barx <- barplot(y.means, names.arg=1:5,ylim=c(0,1.5), col="blue", axis.lty=1, xlab="Replicates", ylab="Value (arbitrary units)")
error.bar(barx,y.means, 1.96*y.sd/10)

y1 <- rnorm(500, mean=1.1)
y1 <- matrix(y1,100,5)
y1.means <- apply(y1,2,mean)
y1.sd <- apply(y1,2,sd)

yy <- matrix(c(y.means,y1.means),2,5,byrow=TRUE)
ee <- matrix(c(y.sd,y1.sd),2,5,byrow=TRUE)*1.96/10
barx <- barplot(yy, beside=TRUE,col=c("blue","magenta"), ylim=c(0,1.5), names.arg=1:5, axis.lty=1, xlab="Replicates", ylab="Value (arbitrary units)")
error.bar(barx,yy,ee)





phdfig22=read.csv(file.choose(), header=T)
phdfig22d<-phdfig22
phdfig22d_long = with(phdfig22d, 
                      rbind( data.frame( lesion=group, stimulus="CS+", latency=CSpluslatency, sem=CSplusSEM ),
                             data.frame( lesion=group, stimulus="CS-", latency=CSminuslatency, sem=CSminusSEM ) )
)
phdfig22d_long$lesion = factor(phdfig22d_long$lesion, levels=unique(phdfig22d_long$lesion))

library(ggplot2)
f3 = ggplot(data = phdfig22d_long, aes(x=lesion, y=latency, fill=stimulus) ) +
  geom_errorbar(aes(ymin=latency, ymax=latency+sem, width = 0.2), position=position_dodge(width=0.90)) + # see http://had.co.nz/ggplot2/position_dodge.html
  geom_bar(stat="identity", position="dodge") + # WORKAROUND: once with no colour, for the legend...
  geom_bar(stat="identity", position="dodge", colour="black", show_guide=FALSE) + # once with a coloured (black!) border, but no legend, or it puts a diagonal line through the legend squares
  
  # see http://groups.google.com/group/ggplot2/browse_thread/thread/bc99c5f7f1c1ab1e/c0b55fa63d42e4e1
  scale_fill_manual(values=c("grey80", "white")) + # grey80 is closer to white than black
  xlab("") +
  scale_y_continuous("Latency to approach (s)", expand=c(0,0), limits = c(0, 5.5) ) +
  theme_bw() 
f3+labs(title = "Autoshaping approach latencies",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(c("left","bottom")), # RNC homebrew to get left/bottom border only
        axis.text.y=element_text(angle=90, hjust=0.5), # rotate the Y axis text anticlockwise 90 degrees, and centre it (0 left, 0.5 centre, 1 right)
        legend.position = c(0.25,0.85),
        legend.key = element_rect() # provides a border to the coloured squares in the legend
)

  