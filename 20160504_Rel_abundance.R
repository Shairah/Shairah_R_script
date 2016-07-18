#read data file rep seq and tax name
getwd
setwd("/Users/shairahar/Dropbox/MiSeqMothur")
repseq <- read.csv("~/Dropbox/MiSeqMothur/repseq.csv") ##either or
repseq=read.csv(file.choose(), header=T)


#make sure the repseq has been modified to resemble the desired number of otu
#the last row should contain the total number of seq of all other otus
repseq <-repseq[-71,]
repseq[70:75,1:5]
repseq[200:205,1:5]
repseq <-repseq[-201,]
repseq <-as.matrix(repseq[,-1], header=F)
dim(repseq) #make sure 3376 otu / 69 samples
head(repseq)
tail(repseq)

TR<-t(repseq) #transpose
rownames(TR) <- paste0("Sample", 1:nrow(TR))
TR[,1:10]

#do group separation
S_T1 = as.matrix(TR[1:6,]) #ones6
S_T2=as.matrix(TR[7:12,]) #ones6
S_T3=as.matrix(TR[13:18,]) #ones6
Sp_T1=as.matrix(TR[19:23,]) #ones5
Sp_T2=as.matrix(TR[24:29,]) #ones6
Sp_T3=as.matrix(TR[30:35,]) #ones6
GW_T1=as.matrix(TR[36:41,]) #ones6
GW_T2=as.matrix(TR[42:45,]) #ones4
GW_T3=as.matrix(TR[46:51,]) #ones6
GWp_T1=as.matrix(TR[52:57,]) #ones6
GWp_T2=as.matrix(TR[58:63,]) #ones6
GWp_T3=as.matrix(TR[64:69,]) #ones6

#create mean vector(1xn)
#need row vector 1 (1xn)
ones6=matrix(1,1,6)
ones5=matrix(1,1,5)
ones4=matrix(1,1,4)
#mean vector
xbST1=(ones6%*%S_T1)*(1/6)
xbST2=(ones6%*%S_T2)*(1/6)
xbST3=(ones6%*%S_T3)*(1/6)
xbSpT1=(ones5%*%Sp_T1)*(1/5)
xbSpT2=(ones6%*%Sp_T2)*(1/6)
xbSpT3=(ones6%*%Sp_T3)*(1/6)
xbGwT1=(ones6%*%GW_T1)*(1/6)
xbGwT2=(ones4%*%GW_T2)*(1/4)
xbGwT3=(ones6%*%GW_T3)*(1/6)
xbGwpT1=(ones6%*%GWp_T1)*(1/6)
xbGwpT2=(ones6%*%GWp_T2)*(1/6)
xbGwpT3=(ones6%*%GWp_T3)*(1/6)

xb.S=rbind(xbST1,xbST2,xbST3)
xb.Sp=rbind(xbSpT1,xbSpT2,xbSpT3)
xb.Gw=rbind(xbGwT1,xbGwT2,xbGwT3)
xb.Gwp=rbind(xbGwpT1,xbGwpT2,xbGwpT3)

#create Diet matrix

Diet=rbind(xb.S,xb.Sp,xb.Gw,xb.Gwp)
Diet[,1:10]
#need to create pivot table
TDiet=t(Diet)
TDiet[1:10,]
View(TDiet[3001:3376,])
# do stuff on X?





#take the first 200 otus
#row 202 and above are composite/total num of seq from other taxa
repseq2 <-as.matrix(repseq[1:201,-1], header=F) 
str(repseq2)
dim(repseq2)
View(repseq2)
head(repseq2)
 
##for plotting## take first 70 otus
repseq3 <-as.matrix(repseq[1:71,-1], header=F) #row 72 and above are composite/toatl num of seq from other taxa
View(repseq3)
dim(repseq3)
rownames(repseq3) <- paste0("OTU", 1:nrow(repseq3))
colnames(repseq3) <- paste0("Sample", 1:ncol(repseq3))

##change data from #seq to frequencies##
freq3=repseq3
for(i in 1:69){freq3[,i]=repseq3[,i]/sum(repseq3[,i])}
View(freq3)

###do nmds and clustering
otutab3 <- t(freq3) 
#calculate distance based on bray-curtis
otutab3.dist=vegdist((otutab3), method='bray', diag= FALSE, upper=FALSE)

#do NMDS ordination
otutab3.mds <- metaMDS(otutab3.dist)
plot(otutab3.dist, type="t") ################?????how to color the graph? or put any character like dot, square etc?
###need to provide label file?
#BS_diet_label

#clustering
#i'd like to see if samples from the same treatment will be clustered together
#k.means clustering
otutab3.kclust=kmeans(otutab3.dist,3,1000)
otutab3.kclust$cluster
View(otutab3.kclust$cluster)




taxa <-read.csv("TaxName.csv", header=T)
head(taxa)
taxmat2 <-as.matrix(taxa[1:71,2:7])
head(taxmat2)
dim(taxmat2)
View(taxmat2)
rownames(taxmat2) <- rownames(repseq3)
rownames(taxmat2) <- rownames(freq3)
class(repseq3)
class(freq3)
class(taxmat2)

##merge both otu and taxonomy table##
#require phyloseq package
library("phyloseq")
REP2 = otu_table(repseq3, taxa_are_rows = TRUE) ##for 200 otus##
REPfreq3=otu_table(freq3, taxa_are_rows = TRUE) ##for 70 otus##
#otu_table works on any numeric matrix. You must also specify if the species are rows or columns
TAX2 = tax_table(taxmat2)
#Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
View(REP2)
str(REP2)
str(TAX2)
View(TAX2)
physeq_gr= phyloseq(REP2,TAX2)
physeq_gr

physeq_fr=phyloseq(REPfreq3,TAX2)
physeq_fr

###create bar plot###
library(ggplot2)
rel.abun1 <- plot_bar(physeq, fill = "Phylum")
rel.abun1 + ggtitle("Relative abundance of Phylum for each sample")

rel.abun2 <- plot_bar(physeq_fr, fill = "Phylum")
rel.abun2 + ggtitle("Relative abundance of Phylum for each sample")

#read Mapping file/labelling#
BLmetadata<-read.csv("~/Dropbox/MiSeqMothur/label.csv")
BLmetadata #The rownames must match the sample names in the otu_table (freq3) if you plan to combine them as a phyloseq-object
freq3
OTUs=freq3
rownames(BLmetadata) = sample_names(physeq_fr)
colnames (OTUs) = sample_names(physeq_fr)
  
OTUs_T1 <- subset(BLmetadata, Time =="1")
OTUs_T2 <- subset(BLmetadata, Time =="2")
OTUs_T3 <- subset(BLmetadata, Time =="3")

Time1 <- OTUs[,BLmetadata[,"Time"]=="1"]
Time2 <- OTUs[,BLmetadata[,"Time"]=="2"]
Time3 <- OTUs[,BLmetadata[,"Time"]=="3"]

S <- OTUs[,BLmetadata[,"Treatment"]=="S"]
Sp <-OTUs[,BLmetadata[,"Treatment"]=="Sp"]
GW <-OTUs[,BLmetadata[,"Treatment"]=="GW"]
GWp <- OTUs[,BLmetadata[,"Treatment"]=="GWp"]



##merge_phyloseq
merge_phyloseq(physeq_fr, BLmetadata)

phyloseq(REP2, TAX2, BLmetadata)

#########
BLrepseq<-as.numeric(BLrepseq)
BLrepseq2<-as.matrix(BLrepseq)
dim(BLrepseq2)
repseq<-as.numeric(BLrepseq [3:3378,2:70])
repseq[1:10, 1:20]
repseq <-as.matrix(repseq)
repseq <-as.numeric(repseq)
repseq <-matrix(repseq, 3376,69)
repseq
str(repseq)
##########

rownames(repseq2) <- paste0("OTU", 1:nrow(repseq2))
colnames(repseq2) <- paste0("Sample", 1:ncol(repseq2))
View(repseq2)
print(repseq2[1:10, 1:5])
      , quote=FALSE)
dim(repseq2)
taxa <-read.csv("TaxName.csv", header=T)
head(taxa)
taxmat <-as.matrix(taxa[1:201,2:7])
head(taxmat)
dim(taxmat)
View(taxmat)
rownames(taxmat) <- rownames(repseq2)
class(repseq2)
class(taxmat)
View(taxmat)
#finish create/read both otu matrix and taxonomymatrix
# In the previous lines, we didn't even need to have phyloseq loaded yet.
# Now we do.
library("phyloseq")
REP = otu_table(repseq2, taxa_are_rows = TRUE) 
#otu_table works on any numeric matrix. You must also specify if the species are rows or columns
TAX = tax_table(taxmat)
#Works on any character matrix. The rownames must match the OTU names (taxa_names) of the otu_table if you plan to combine it with a phyloseq-object.
View(REP)
str(REP)
str(TAX)
View(TAX)
physeq = phyloseq(REP,TAX)
physeq

#plotting relative abundance
library(ggplot2)
data("sampledata")
plot_bar(physeq1, fill = "Phylum")




##metadata
sampledata <- read.csv("~/Dropbox/MiSeqMothur/label.csv", header=T)
View(sampledata)
rownames(sampledata) <- paste0("Sample", 1:ncol(repseq2))
row.names(sampledata) <-  sample_names(physeq)

library("ape")
random_tree = rtree(69, rooted = TRUE, tip.label = sample_names(physeq))
plot(random_tree)
 #69 is the number of samples
physeq1 = merge_phyloseq(physeq, sampledata)
physeq1
plot_heatmap(physeq1)
plot_heatmap(physeq1, taxa.label = "Phylum")


physeq2 = phyloseq(REP, TAX, sampledata, random_tree)


# phyloseq-class experiment-level object
## otu_table()   REP Table:         [ 10 taxa and 10 samples ]
## sample_data() Sample Data:       [ 10 samples by 2 sample variables ]
## tax_table()   Taxonomy Table:    [ 10 taxa by 7 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 10 tips and 9 internal nodes ]


dev.off()

####Time 1####
relabund1T1dup <- read.csv("~/Dropbox/MiSeqMothur/relabund1T1-2.csv", header=T)
relabund1T1<-relabund1T1dup[,-1]
relabund1T1<-relabund1T1[,-1]
relabund1T1dup[,1]
relabund1T1<-as.matrix(relabund1T1)
str(relabund1T1)
Comm=rbind("S","Sp", "Gw", "Gwp", "Stream", "GWTR")
Comm
colors<-c("#b3cde3","black","#fbb4ae","black","white","black"),lwd=1, angle=45, density=c(0,5), add=T))

barplot(relabund1T1,beside=T, col=c("#b3cde3","#b3cde3","#fbb4ae","#fbb4ae","white","white"), 
  ylim=c(0,0.3), ylab=("Relative abundance"), xaxt="n", xlab=" ", main="Relative abundance of taxa for fish gut and water communities at Time 1")
barplot(relabund1T1,beside=T, col=c("#b3cde3","black","#fbb4ae","black",    "white","black"), 
        legend= rownames(relabund1T1), yaxt="n", xlab=" ", las=2, cex.names=0.65, lwd=1, angle=45, density=c(0,20,0,20,0,20), add=T)

counts <- table(mtcars$gear)
barplot(counts, main="Car Distribution by Gears and VS",
        xlab="Number of Gears", col=c("#7fc97f","#beaed4"),
        legend = rownames(counts),beside=T)
barplot(counts, col=c("black"),
        legend = rownames(counts), beside=T, lwd=1, angle=45, density=c(0,5), add=T)

dev.off()

####Time 3 0.1% taxa####
relabund0.1T3 <- read.csv("~/Dropbox/MiSeqMothur/relabund0.1t3.csv", header=T)
head(relabund0.1T3)
relabund0.1T3new <- relabund0.1T3[,c("Trt","Genera","prelab0.1t3","psem")]
head(relabund0.1T3new)
#subset for water, food, gut #https://dzone.com/articles/learn-r-how-create-data-frames
foodT3<-subset(relabund0.1T3new[1:32,])
str(foodT3)
waterT3 <- subset(relabund0.1T3new[33:64,])
gutT3<-subset(relabund0.1T3new[65:128,])
str(waterT3) 
#drop level that do not occur
waterT3$Trt <- factor(waterT3$Trt)
waterT3$Genera <- factor(waterT3$Genera)
foodT3$Trt <- factor(foodT3$Trt)
foodT3$Genera <- factor(foodT3$Genera)
gutT3$Trt <- factor(gutT3$Trt)
gutT3$Genera <- factor(gutT3$Genera)


# order things based on the order seen in the data frame:
x$variable <- factor(x$variable, levels=unique(as.character(x$variable)) )
#order levels based on another variable (value in this case):
x <- transform(x, variable=reorder(variable, -value) ) 

#plot water
library(ggplot2)
#reorder to make sure genus appear as sorted
waterT3$Genera <- factor(waterT3$Genera, levels=unique(as.character(waterT3$Genera)) )
# Default bar plot
wT3<- ggplot(waterT3, aes(x=Genera, y=prelab0.1t3, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t3, ymax=prelab0.1t3+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(waterT3$prelab0.1t3), max(waterT3$prelab0.1t3), by = 5)), expand=c(0,0))

# Finished bar plot
wT3+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c('black','white')) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))

#plot food
library(ggplot2)
foodT3$Genera <- factor(foodT3$Genera, levels=unique(as.character(foodT3$Genera)) )
# Default bar plot
fT3<- ggplot(foodT3, aes(x=Genera, y=prelab0.1t3, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t3, ymax=prelab0.1t3+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(foodT3$prelab0.1t3), max(foodT3$prelab0.1t3), by = 5)), expand=c(0,0))

# Finished bar plot
fT3+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))

#plot gut
library(ggplot2)
gutT3$Genera <- factor(gutT3$Genera, levels=unique(as.character(gutT3$Genera)) )
# Default bar plot
gT3<- ggplot(gutT3, aes(x=Genera, y=prelab0.1t3, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t3, ymax=prelab0.1t3+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(gutT3$prelab0.1t3), max(gutT3$prelab0.1t3), by = 5)), expand=c(0,0))

# Finished bar plot
gT3+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c("hotpink","red3","cyan","blue")) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))

####Time 2 0.1% taxa####
relabund0.1T2 <- read.csv("~/Dropbox/MiSeqMothur/relabund0.1t2.csv", header=T)
head(relabund0.1T2)
relabund0.1T2new <- relabund0.1T2[,c("Trt","Genera","prelab0.1t2","psem")]
head(relabund0.1T2new)
#subset for water, food, gut #https://dzone.com/articles/learn-r-how-create-data-frames
foodT2<-subset(relabund0.1T2new[1:32,])
str(foodT2)
waterT2 <- subset(relabund0.1T2new[33:64,])
gutT2<-subset(relabund0.1T2new[65:128,])
str(waterT2) 
#drop level that do not occur
waterT2$Trt <- factor(waterT2$Trt)
waterT2$Genera <- factor(waterT2$Genera)
foodT2$Trt <- factor(foodT2$Trt)
foodT2$Genera <- factor(foodT2$Genera)
gutT2$Trt <- factor(gutT2$Trt)
gutT2$Genera <- factor(gutT2$Genera)
#plot water
library(ggplot2)
#reorder to make sure genus appear as sorted
waterT2$Genera <- factor(waterT2$Genera, levels=unique(as.character(waterT2$Genera)) )
# Default bar plot
wT2<- ggplot(waterT2, aes(x=Genera, y=prelab0.1t2, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t2, ymax=prelab0.1t2+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(waterT2$prelab0.1t2), max(waterT2$prelab0.1t2), by = 5)), expand=c(0,0))

# Finished bar plot
wT2+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c('black','white')) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))

#plot food
library(ggplot2)
foodT2$Genera <- factor(foodT2$Genera, levels=unique(as.character(foodT2$Genera)) )
# Default bar plot
fT2<- ggplot(foodT2, aes(x=Genera, y=prelab0.1t2, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t2, ymax=prelab0.1t2+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(foodT2$prelab0.1t2), max(foodT2$prelab0.1t2), by = 5)), expand=c(0,0))

# Finished bar plot
fT2+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))

#plot gut
library(ggplot2)
gutT2$Genera <- factor(gutT2$Genera, levels=unique(as.character(gutT2$Genera)) )
# Default bar plot
gT2<- ggplot(gutT2, aes(x=Genera, y=prelab0.1t2, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t2, ymax=prelab0.1t2+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(gutT2$prelab0.1t2), max(gutT2$prelab0.1t2), by = 5)), expand=c(0,0))

# Finished bar plot
gT2+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c("hotpink","red3","cyan","blue")) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))

####Time 1 0.1% taxa####
relabund0.1T1 <- read.csv("~/Dropbox/MiSeqMothur/relabund0.1t1.csv", header=T)
head(relabund0.1T1)
relabund0.1T1new <- relabund0.1T1[,c("Trt","Genera","prelab0.1t1","psem")]
head(relabund0.1T1new)
#subset for water, gut #https://dzone.com/articles/learn-r-how-create-data-frames
waterT1 <- subset(relabund0.1T1new[1:32,])
gutT1<-subset(relabund0.1T1new[33:96,])
str(waterT1) 
#drop level that do not occur
waterT1$Trt <- factor(waterT1$Trt)
waterT1$Genera <- factor(waterT1$Genera)

gutT1$Trt <- factor(gutT1$Trt)
gutT1$Genera <- factor(gutT1$Genera)
#plot water
library(ggplot2)
#reorder to make sure genus appear as sorted
waterT1$Genera <- factor(waterT1$Genera, levels=unique(as.character(waterT1$Genera)) )
# Default bar plot
wT1<- ggplot(waterT1, aes(x=Genera, y=prelab0.1t1, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t1, ymax=prelab0.1t1+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(waterT1$prelab0.1t1), max(waterT1$prelab0.1t1), by = 5)), expand=c(0,0))

# Finished bar plot
wT1+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c('black','white')) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))


#plot gut
library(ggplot2)
gutT1$Genera <- factor(gutT1$Genera, levels=unique(as.character(gutT1$Genera)) )
# Default bar plot
gT1<- ggplot(gutT1, aes(x=Genera, y=prelab0.1t1, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=prelab0.1t1, ymax=prelab0.1t1+psem), width=.2,
                position=position_dodge(.9)) + 
  scale_y_continuous(breaks = round(seq(min(gutT1$prelab0.1t1), max(gutT1$prelab0.1t1), by = 5)), expand=c(0,0))

# Finished bar plot
gT1+labs( x="Genera", y = "Relative abundance (%)")+
  theme_classic() +
  scale_fill_manual(values=c("hotpink","red3","cyan","blue")) +
  theme(axis.text.x = element_text(angle = 45, hjust=0.85, size = 12), axis.text.y = element_text(size =12))




#use base package to create bar plot
maxplotTop <- max(gutT3$prelab0.1t3) +
  gutT3[gutT3$prelab0.1t3 == max(gutT3$prelab0.1t3), 6]*3
basegT3 <- barplot(height = gutT3$prelab0.1t3,
                      beside = TRUE, las = 1,
                    
                      cex.names = 0.75,
                      main = "Rel. abundance of predominance genera for gut communities at Time3",
                      ylab = "Relative abundance (%)",
                      xlab = "Genera",
                      border = "black", axes = TRUE,
                      legend.text = TRUE,
                      args.legend = list(title = "Group", 
                                         x = "topright",
                                         cex = .7))

segments(basegT3, gutT3$prelab0.1t3 , basegT3,
         gutT3$prelab0.1t3 + gutT3$psem, lwd = 1.5)

arrows(basegT3, gutT3$prelab0.1t3 , basegT3,
       gutT3$prelab0.1t3 + gutT3$psem, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)

text(x = basegT3, y = par("usr")[3] - 1, srt = 45,
     adj = 1, labels = gutT3$Genera, xpd = TRUE)


#####################

relabund1T1<-relabund1T1dup[,-1]
relabund1T1<-relabund1T1[,-1]
relabund1T1dup[,1]
relabund1T1<-as.matrix(relabund1T1)
str(relabund1T1)
relab0.1t3



text(relabund1T1, labels = xlab, las=2, adj = 1, xpd = TRUE, cex.axis=0.8)

Otuname<-as.vector(colnames(relabund1T1))

par(mar= c(7,4,2,2) +0.2)
end_point = 0.5 + nrow(relabund1T1) + nrow(relabund1T1)-1
text(1:8, par("usr")[3] - 0.25, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
axis(2)
?par
text(seq(1.5,end_point,by=2), par("usr")[3]-0.25, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(rownames(mtcars)), cex=0.65)
text(relabund1T1,par("usr")[3]-0.25, srt=60, adj=1,
     labels=Otuname,xpd=TRUE)
text(seq(1.5,end_point,by=2), par("usr")[3]-0.25, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = Otuname, cex=0.65)
mtext(1, text = "X Axis Label", line = 6)

## Increase bottom margin to make room for rotated labels
par(mar = c(7, 4, 4, 2) + 0.1)
## Create plot with no x axis and no x axis label
plot(1 : 8, xaxt = "n",  xlab = "")
## Set up x axis with tick marks alone
axis(1, labels = FALSE)
## Create some text labels
labels <- paste("Label", 1:8, sep = " ")
## Plot x axis labels at default tick marks
text(1:8, par("usr")[3] - 0.25, srt = 45, adj = 1,
     labels = labels, xpd = TRUE)
## Plot x axis label at line 6 (of 7)
mtext(1, text = "X Axis Label", line = 6)



axis(1, at=seq(1, 10, by=1), labels = FALSE)
text(seq(1, 10, by=1), par("usr")[3] - 0.2, labels = Otuname, srt = 45, pos = 1, xpd = TRUE)




colnames(relabund1T1)
Otuname<-colnames(relabund1T1[,-1])
Otuname


library(ggplot2)
relabund1T1b<- ggplot(relabund1T1, aes(x=Otus, y=relab1, fill=Trt)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9)) 
relabund1T1b+labs(title="Relative abundance of OTU more than 1%", x="Genera", y = "Rel. abund")+
  theme_classic() +
  scale_fill_manual(values=c("#b3cde3","#b3cde3","#fbb4ae","#fbb4ae","white","white")) +
  scale_fill_manual(values=c("#b3cde3","black","#fbb4ae","black","white","black"),  
 lwd=1, angle=45, density=c(0,15,0,15,0,15), add=T)

print(relabund1T1b)
