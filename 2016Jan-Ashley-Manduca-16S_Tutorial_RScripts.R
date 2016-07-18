#Day 3 new commands and scripts
library(vegan)
otu=read.table("Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)

#Before beginning analysis, check that the OTU table looks as we expect:
head(otu)

row.names(otu)
colnames(otu)

dim(otu)
rdp=otu[,77]
names(rdp)=row.names(otu)
otu=otu[,-77]
dim(otu)

#Read in mapping file
map=read.table("Manduca_map_sorted.txt", header=TRUE, sep="\t")
#attach(map)

  
#Make a bray-curtis resemblance matrix
braycurtis.d=vegdist(t(otu), method="bray")

#make a sorenson's resemblance matrix
sorenson.d=vegdist(t(otu), method="bray", binary=TRUE)

#Read in the phylogenetic matrices
weighted_u=read.table("weighted_unifrac_Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t")
#attach(map)
head(weighted_u)
weighted_u.d=as.dist(weighted_u)
head(weighted_u.d)

unweighted_u=read.table("unweighted_unifrac_Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t")
#attach(map)
head(unweighted_u)
unweighted_u.d=as.dist(unweighted_u)
head(unweighted_u.d)

#make a heatmap
library(gplots)
colors <- colorRampPalette(c("white", "black"))(20)
heatmap.2(t(otu), trace="none", col=colors)


#make a quick NMDS ordination
braycurtis.mds <- metaMDS(braycurtis.d)
plot(braycurtis.mds, type="t")

#Perform all nmds and plot at once
sorenson.mds=metaMDS(sorenson.d)
weighted_u.mds=metaMDS(weighted_u.d)
unweighted_u.mds=metaMDS(unweighted_u.d)

#plot all the ordinations to compare
par(mfrow=c(2,2))
plot(braycurtis.mds, type="t", main= "Bray-Curtis")
plot(sorenson.mds, type="t", main="Sorenson")
plot(weighted_u.mds, type="t", main="Weighted UniFrac")
plot(unweighted_u.mds, type="t", main="Unweighted UniFrac")


#write out tables
braycurtis=as.matrix(braycurtis.d)
write.table(braycurtis, "BrayCurtis_even861.txt",row.names=TRUE, col.names=TRUE,sep="\t", quote=FALSE)

sorenson=as.matrix(sorenson.d)
write.table(sorenson, "Sorenson_even861.txt",row.names=TRUE, col.names=TRUE,sep="\t", quote=FALSE)

#Day 4 New Scripts
#Read in tables:
unweighted_u=read.table("unweighted_unifrac_Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t")
unweighted_u.d=as.dist(unweighted_u)

weighted_u=read.table("weighted_unifrac_Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t")
weighted_u.d=as.dist(weighted_u)

braycurtis=read.table("BrayCurtis_even861.txt", header=TRUE, row.names=1, sep="\t")
braycurtis.d=as.dist(braycurtis)

sorenson=read.table("Sorenson_even861.txt", header=TRUE, row.names=1, sep="\t")
sorenson.d=as.dist(sorenson)

otu=read.table("Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
#Don't forget to remove rdp ids:
rdp=otu[,77]
names(rdp)=row.names(otu)
otu=otu[,-77]

map=read.table("Manduca_map_sorted.txt", header=TRUE, sep="\t")

#Load libraries
library(vegan)
library(limma)

#Hypothesis testing.
#permutated analysis of variance (PERMANOVA)
#d is the input distance matrix, and A is the grouping of interest
ad=adonis(braycurtis.d~Treatment, data=map, permutations=999)
a.table=ad$aov.tab

#Permutated multivariate analysis of beta-dispersion
b=betadisper(braycurtis.d, group=map[,"Treatment"], type="median")
b
b.perm=permutest(b, group=map[,"Treatment"], type="median", permutations=999, pairwise=TRUE)
b.perm

#How reproducible are replicates?
#1.  Make a vector of all of the replicate names
u=unique(map[,"Replicates"])

#Make a loop to cycle through each set of replicates.  
meanreps.out=NULL
for(i in 1:length(u)){
  
  #first, reduce the resemblance matrix (table version) to only the data you need for each rep
  temp=braycurtis[map[,"Replicates"]==u[i],map[,"Replicates"]==u[i]]
  
  #then, convert thet able to an R-recognized resemblance vector using as.dist()
  temp.d=as.dist(temp)
  
  #Calculate the mean of the vector
  m=mean(temp.d)
  
  #Add it to a growing vector of means
  meanreps.out=c(meanreps.out,m)
}

#Apply the names to the vector, from u
names(meanreps.out)=u

#Plot the results
hist(meanreps.out)

##Make a time distance matrix
makeTimeDist.f=function(map_file){
  map=map_file
  temp=as.matrix(map[,"Instar"])
  names(temp)=map[,"SampleID"]
  temp.d=dist(temp, method="manhattan", diag=FALSE)
  head(temp.d)
  temp.out=as.matrix(temp.d)
  colnames(temp.out)=map[,"SampleID"]
  row.names(temp.out)=map[,"SampleID"]
  print(head(temp.out))
  return(temp.d)
}

time.d=makeTimeDist.f(map)

mantel(time.d,braycurtis.d, method="pearson", permutations = 999)
mantel(time.d,sorenson.d, method="pearson", permutations = 999)
mantel(time.d,weighted_u.d, method="pearson", permutations = 999)
mantel(time.d,unweighted_u.d, method="pearson", permutations = 999)

#Fitting environmental vectors to axis scores, and an example of Correspondance analysis
otu.ca=cca(t(otu))
plot(otu.ca)

#Use envfit to fit the instar information to the CA.
ev.instar=envfit(otu.ca, as.numeric(map[,"Instar"]), perm=1000)
ev.instar
plot(otu.ca,display="site")
plot(ev.instar)

#using PERMANOVA to understand time-treatment interactions:
ad2=adonis(braycurtis.d~Treatment*as.numeric(Instar), data=map, permutations=999)
a.table2=ad2$aov.tab
a.table2



###New Scripts Day 5.  Use previous scripts to read in all needed files.
## Read in files:
braycurtis=read.table("BrayCurtis_even861.txt", header=TRUE, row.names=1, sep="\t")
braycurtis.d=as.dist(braycurtis)
otu=read.table("Manduca_otu_table_even861_sorted.txt", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
rdp=otu[,77]
names(rdp)=row.names(otu)
otu=otu[,-77]

map=read.table("Manduca_map_sorted.txt", header=TRUE, sep="\t")

#PROTEST
#Designate treatment groups using the unique() function
u=unique(map[,"Treatment"])

#Make a resemblence matrix for each treatment group
makeRedresem.f=function(resem_fp,group){
  resem=resem_fp
  red=resem[map[,"Treatment"]==group, map[,"Treatment"]==group]
  red.d=as.dist(red)
  return(red.d)
}

#note that we use the "table" form of the resemblance matrix, not the vector for, as input for this function
#Also note that this only works if there are the EXACT SAME number of samples in the EXACT SAME order.
mrs.d=makeRedresem.f(braycurtis,"MRS")
pbs.d=makeRedresem.f(braycurtis,"PBS")
LGG.d=makeRedresem.f(braycurtis,"LGG")
Lp.d=makeRedresem.f(braycurtis,"Lp")
Lr.d=makeRedresem.f(braycurtis,"Lr")

#Run each pair test between all pairs of treatments, e.g.:
protest(mrs.d,LGG.d, permutations=999)
protest(Lr.d,LGG.d, permutations=999)

#Also with mantel
mantel(mrs.d, LGG.d, method="pearson", permutations=999)
mantel(Lr.d, LGG.d, method="pearson", permutations=999)

#Venn analysis
#Make a concatenated dataset by treatment, and calculate shared OTUs (taxonomic Venn analysis)
#concatenate dataset so that all of the sequences from the same day are together
otu.trt=matrix(0,ncol=length(u),nrow=nrow(otu))
for(i in 1:length(u)){
  temp=otu[,map[,"Treatment"]==u[i]]
  #print(i)
  dim(temp)
  temp2=rowSums(temp)
  otu.trt[,i]=temp2
}
row.names(otu.trt)=row.names(otu)
colnames(otu.trt)=u

#Make VennCounts

#Install limma
#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")

library(limma)
v=vennCounts(otu.trt)
#Write out the results of venncounts
write.table(v, "VennCounts.txt", quote=FALSE, sep="\t")

#What is the taxonomic identity of  the OTUs in each Venn category?
#Make a presence/absence (binary) table
otu.trt.pa=1*(otu.trt>0)
#u=unique(map[,"Treatment"])

#This loop will output one text file per VennCounts category that has a list of all of the OTU IDs and their RDP taxonomic affiliations that belong to each VennCount group.  
#The group name is the row of the Venn Count, so you will have to compare the VennCount file to the output to determine which group.
tmp=NULL
#VennOTUGroups=NULL
VennOTUGroups=vector(length=nrow(otu))

for(y in 1:nrow(v)){
  if(v[y,ncol(v)]!=0){
    
    for(z in 1:nrow(otu.trt.pa)){
      
      if(sum(1*(as.vector(otu.trt.pa[z,]) == as.vector(v[y,(1:ncol(v)-1)])))==length(u)){
        n=c(row.names(otu.trt.pa)[z],paste(rdp[z]))
        tmp=rbind(tmp,n)
        #VennOTUGroups=c(VennOTUGroups,y)
        VennOTUGroups[z]=y
      }

    }
    
    write.table(tmp, paste("VennCounts_Condition_",y,"_OTUList.txt",sep=""), row.names=FALSE, sep="\t", quote=FALSE)
    print(dim(tmp))
  }
  tmp=NULL
}

#Check VennOTUGroups.  This can be added to the mapping file and used for future exploration.
#It should be the match the number of OTUs (no. rows in your otu table), and should include numbers from 1 to the number of VennGroups
#In this case, with six starting treatment groups, we have 64 VennGroups of overlapping taxa between all combinations
length(VennOTUGroups)

#Check that VennOTUGroups contains as you expect from v, the VennCounts result.
#Using this logical statement, we sum the number of times that it is TRUE that VennOTUGroups 
#has a 2.  From the v, we know that there should be 63 OTUs that are only viewed in MRS, which is VennOTUGroups 2.
#We have a match!
#Check on this?
sum(1*(VennOTUGroups==2))
head(VennOTUGroups)

#Clustering OTUs by similar occurence patterns
#Make a table relativized by rows:
rSums=rowSums(otu)
otu.relrows=otu/rSums

#use the hclust function to apply hierarchical clustering
braycurtis.rows.d=vegdist(otu.relrows, method="bray")
otu.cluster=hclust(braycurtis.rows.d, method="complete")
plot(otu.cluster)

#We also can use our VennOTUGroups to find if there are clusters of OTUs with shared presence/absences
plot(otu.cluster,labels=VennOTUGroups)

#We use make.cephnames() to shorten the rdp
rdp2=make.cepnames(rdp)
plot(otu.cluster, labels=rdp2)

#Making SADs
#Model fitting to species abundance to the rowSums of the OTU table, inspect R.  You want to chose the model with the lowest AIC and BIC
#(In actuality, none of these models are spectacular, but Mandelbrot is best)
r.sad=radfit(rSums)
r.sad
#here, it is clear from the plot that the Mandelbrot model fits best
quartz()
plot(r.sad, las=1,main = "Manduca species abundance distribution", ylab="Abundance (No. sequences)", xlab="OTUs ranked by abundance")

#Species occurrence distributions- based on the binary matrix
#make a binary matrix
otu.pa=1*(otu>0)
rSums.pa=rowSums(otu.pa)

#Here, we see that a few models are equally acceptable
r.sod=radfit(rSums.pa)
r.sod
plot(r.sod, las=1,ylab="Occurrence (out of 76 total observations)", xlab="OTUs ranked by occurrence", main="Manduca species occurrence distribution")


#MultiCOLA
##Mantel MultiCOLA with Bray-Curtis distances
MantelMultiCOLA.f=function(otu_fp){
  #Step 1.  Read in full dataset:
  otu2=otu_fp
  
  library(vegan)
  
  cutoff=c(1.00,0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05, 0.025, 0.01, 0.005, 0.001)
  
  m.out=NULL
  
  otu.pa=1*(otu2>0)
  
  r=rowSums(otu2)
  
  #Create a vector of indexes of the ranked OTUs
  r2=sort(r,decreasing=TRUE, index=TRUE)
  r2$ix
  r2$x
  
  for(j in 1:length(cutoff)){
    print(cutoff[j])
    
    no.keep=ceiling(cutoff[j]*nrow(otu2))
    otu.index.keep=r2$ix[1:no.keep]
    otu.keep=otu2[otu.index.keep,]
    print(head(otu.keep))
    
    
    #Write out otu tables at each cutoff
    #write.table(otu.keep, paste(u[i],"_",cutoff[j],"_otu.txt", sep=""),sep="\t", quote=FALSE)
    
    all.dist=vegdist(t(otu2), method="bray")
    subset.dist=vegdist(t(otu.keep),method="bray")
    
    m1=mantel(all.dist,subset.dist, method="pearson",permutations=999)
    
    
    m=c(paste(cutoff[j]),m1$statistic,m1$signif,dim(otu.keep)[1])
    m.out=rbind(m.out,m)
  }  	
  
  colnames(m.out)=c("Cutoff", "AllvSubsetPearsonR", "AllvSubset_pvalue", "NoOTUsSubset")
  
  #write.table(m.out, "MantelMultiCOLA.txt", sep="\t", quote=FALSE, row.names=FALSE)
  return(m.out)
}

#Run the function:
MultiCOLA.test=MantelMultiCOLA.f(otu)

#Inspect the output.  From this, we see that there is one OTU that is driving the dynamics that we see.
#We can deduce this because the most reduced dataset, including only 1 taxon comprising the upper 0.001 percent abundance of all taxa, is still highly correlated (0.88) to the original
MultiCOLA.test

#Who is the lone taxon that is so powerful??  Check our rdp
rdp["248"]

otu248=as.numeric(otu["248",])


####End Day 5 functions

########Misc useful R commands and functions ###############
#Relativize the OTU table
#make relative abundance table
makeRAtable.f=function(data){
  cSum1<-colSums(data)
  
  #define an empty matrix to put your RF values into
  newdata<-matrix(0,dim(data)[1], dim(data)[2])
  
  #Assign the same column and row names to the new matrix.
  colnames(newdata)<-colnames(data)
  rownames(newdata)<-rownames(data)
  
  #Each cell will be divided by the column sum.
  for (i in 1:length(data)){
    newdata[,i] <- data[,i]/cSum1[i]
  }
  
  return(newdata)
}

#Removing singletons from the dataset:
otu.nosigs=otu[rowSums(otu)>1,]

#Removing doubletons from the dataset:
otu.nodoubs=otu[rowSums(otu)>2,]

#Other hypothesis tests for differences between groups

#multi-response permutation procedure (MRPP)
m=mrpp(braycurtis.d,grouping=map[,"Treatment"], permutations=999)
m

#analysis of similarity (ANOSIM)
as=anosim(braycurtis.d, grouping=map[,"Treatment"], permutations=999)
as

otu248=as.numeric(otu["248",])
otu248.LGG=otu248[map[,"Treatment"]=="LGG"]
otu248.Lr=otu248[map[,"Treatment"]=="Lr"]
otu248.Lp=otu248[map[,"Treatment"]=="Lp"]
otu248.PBS=otu248[map[,"Treatment"]=="PBS"]

test=aov(otu248~Treatment, data=map)

test2=cca(X=t(otu), Z=otu248)