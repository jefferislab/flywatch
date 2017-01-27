#This is code to extract and plot the the Distances from the center
library(tiff)
library(rmngb)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xlsx)
library(PMCMR)
library(car)

#Load up required functions
Q1Q4.DistanceCenter<-function(data) {
  Q1Q4.DistanceCenter<-matrix(ncol=2, byrow=TRUE, data=c(1,20))
  output<-sapply(data, "[", Q1Q4.DistanceCenter)
  output
}
Q2Q3.DistanceCenter<-function(data) {
  Q2Q3.DistanceCenter<-matrix(ncol=2, byrow=TRUE, data=c(1,22))
  output<-sapply(data, "[",  Q2Q3.DistanceCenter)
  output
}
minimum.frames<-function(x) min(apply(X = x,1, length)) #Determine the final frame as there is variation
frame.cut<-function(mat) mat[,1:finalframe]
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
minimum.frames<-function(x) min(apply(X = x,1, length)) #Determine the final frame as there is variation
frame.cut<-function(mat) mat[,1:finalframe]

#Sort through the directory and pull out the TIFFs.
dir<-list.files(pattern="measurementsXY.tif$", recursive=TRUE) #Pull out the correct tif
#Find the different genotypes by pulling out all the genotype data from TIFF filenames
if(length(unique(sapply(strsplit(dir, "_"), "[", 7)))==1) {
  genotypes<-sapply(strsplit(dir, "_"), "[", 6)
  ugenotypes<-unique(genotypes)
  ugenotypes<-ugenotypes[(-1)*c(2, 5:8, 18, 19, 33,46,50)]
}
if(length(unique(sapply(strsplit(dir, "_"), "[", 7)))>1) {
  genotypes<-sapply(strsplit(dir, "_"), "[", 7)
  genotypes<-sapply(strsplit(genotypes, "/"), "[", 1)
  ugenotypes<-unique(genotypes)[-1]
}
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes)

#Initialise and run the data extraction.
PI<-vector("list", length = length(ugenotypes));names(PI)<-ugenotypes
Q1Q4.Dist<-vector("list", length = length(ugenotypes));names(Q1Q4.Dist)<-ugenotypes
Q2Q3.Dist<-vector("list", length = length(ugenotypes));names(Q2Q3.Dist)<-ugenotypes
for(i in 1:length(ugenotypes)) {
  files<-grep(paste0("_", ugenotypes[i]),
              dir, value=TRUE, fixed=TRUE)
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    Q1Q4.Dist[[i]]<-rbind(Q1Q4.Dist[[i]],Q1Q4.DistanceCenter(data))
    Q2Q3.Dist[[i]]<-rbind(Q2Q3.Dist[[i]],Q2Q3.DistanceCenter(data))
  }
}


#Fix the length of the frames for further analysis
finalframe<-min(melt(lapply(Q1Q4.Dist, FUN = minimum.frames))[,1])
Q1Q4.Dist<-lapply(Q1Q4.Dist, frame.cut)

finalframe<-min(melt(lapply(Q2Q3.Dist, FUN = minimum.frames))[,1])
Q2Q3.Dist<-lapply(Q2Q3.Dist, frame.cut)

index<-c(1:finalframe)/30 #This is the framerate from the camerasettings.json file

#Dist is calculated as (D-D0)/D0. Note that Q3Q4 stimulation happens first!
single_dist<-function(x) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(x)))
  First30<-colMeans(m[m$seconds<=30 ,])[-1]
  Stim1<-colMeans(m[m$seconds>=30 & m$seconds<=60 ,])[-1]
  Stim2<-colMeans(m[m$seconds>=30 & m$seconds<=60 ,])[-1]
  data.frame(First30, Stim1, Stim2)
}
mat_single_dist<-function(x) apply(x, 1, single_dist) #These function calculate average dist proxy for each genotype at stim and control timepoints
Q1Q4single<-melt(sapply(Q1Q4.Dist, mat_single_dist)) #All the averages for Q1Q4 dist at control and stim timepoints
Q1Q4single$quad<-"Q1Q4"
Q2Q3single<-melt(sapply(Q2Q3.Dist, mat_single_dist)) #All the averages for Q2Q3 dist at control and stim timepoints
Q2Q3single$quad<-"Q2Q3"


dist<-rbind(Q1Q4single, Q2Q3single)
names(dist)<-c("Stim", "value", "N", "Genotype", "quad")
dist_singleVal<-data.frame()
for(i in 1:length(ugenotypes)) {
  data<-filter(dist, Genotype==ugenotypes[i]) #First select everything with the same genotype

  #calculate D0, all quads and the First30
  D0<-filter(data, Stim=="First30")
  D0<-aggregate(x = D0$value, by=list( D0$N), FUN=mean)
  names(D0)<-c("N", "value")

  #caluclate D, which is the average of the stimulation quadrants' averages
  D<-rbind(filter(data, Stim=="Stim1" & quad=="Q2Q3"), filter(data,Stim=="Stim2" & quad=="Q1Q4"))
  D<-aggregate(x = D$value, by=list( D$N), FUN=mean)
  names(D)<-c("N", "value")

  #Perform the delta L calculation
  deltaD<-((D-D0)/D0)[,-1] #This is a df of all the single value differences in locomotion proxy. One for each experiment.
  deltaD<-data.frame(Genotype=ugenotypes[i], deltaD)
  dist_singleVal<-rbind(dist_singleVal, deltaD)
}
dist_singleVal<-arrange(dist_singleVal, desc(Genotype=="EmptySp"))
dist_singleVal[,1]<-as.factor(as.character(dist_singleVal[,1])) #Reorganise the factor levels as dunn test function takes first level
save(dist_singleVal, file=paste0(getwd(),"/dist_data.rda"))

#Determine significance and plot the data
#Function to pull out the p-values
pvals<-as.data.frame(dunn.test.control(x = dist_singleVal$deltaD,
                                       g= as.factor(dist_singleVal$Genotype), p.adjust.method = "fdr")$p.value)
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  dist_singleVal, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.1, "Significant", "Not Significant")


g<-ggplot(data=pvals, aes(x=reorder(Genotype, deltaD), y=deltaD))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="deltaD",title="") #Titles
g<-g+theme(legend.title=element_blank())
g


#After running the analyses for all the screen data, let's put everything together for distance
#ANALYSIS III: Combining data from both screening sessions, plotting both and analysing repeats
setwd("/Volumes/Data/BehaviourData/")
Dec2015<-loadRData("./Mike_newrig_Dec2015_screen/dist_data.rda")
Dec2015<-filter(Dec2015, Genotype!="Empty") #Bug in the DunnTest code, need to remove this to use EmptySp as control
Dec2015<-filter(Dec2015, Genotype!="Gr66a")
Sept2016<-loadRData("./Mike_newrig_Sept2016_screen/dist_data.rda")
Sept2016<-filter(Sept2016, Genotype!="Gr66a")
all.data<-rbind(Dec2015,Sept2016)
#Remove unwanted samples
all.data<-filter(all.data, Genotype!="test"| Genotype!="MB83c")
all.data<-arrange(all.data, desc(Genotype=="EmptySp"))
all.data[,1]<-as.factor(as.character(all.data[,1])) #Reorganise the factor levels as dunn test function takes first level

pvals<-as.data.frame(dunn.test.control(x = all.data$deltaD,
                                       g= as.factor(all.data$Genotype), p.adjust.method = "fdr")$p.value)
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  all.data, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.1, "Significant", "Not Significant")

g<-ggplot(data=pvals, aes(x=reorder(Genotype, deltaD, median), y=deltaD))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="deltaL",title="") #Titles
g<-g+theme(legend.title=element_blank())
g
ggsave(filename = "allscreen_distance_from_center_proxy_FDR10.pdf", width=16
       ,height =9 ,plot=g, path=".")

