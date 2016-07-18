
map[,"Classification"]
Logic<-map[,"Classification"]=="Recovered"
Logic
Recovered_OTU_Table=soils[,Logic]

Recovered_OTU_Table=soils[,map[,"Classification"]=="Recovered"]


map<-read.table("FishMetadata.txt")
OTUs<- read.table("SA.txt")
OTUs_T1 <- OTUs[,map[,"Treatment"]=="1"]