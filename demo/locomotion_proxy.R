#This is code to extract and plot the locomoter proxy
library(tiff)
library(rmngb)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xlsx)
library(PMCMR)
library(car)

#Load up required functions
PIextract<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2)) #Creates a matrix to give position for the PI value data using subset, row1 col2
  PItime<-sapply(data, "[", s)
  PItime
}
Q1Q4.locomotion<-function(data) {
  Q1Q4.locomotion<-matrix(ncol=2, byrow=TRUE, data=c(1,16))
  output<-sapply(data, "[", Q1Q4.locomotion)
  output
}
Q2Q3.locomotion<-function(data) {
  Q2Q3.locomotion<-matrix(ncol=2, byrow=TRUE, data=c(1,17))
  output<-sapply(data, "[", Q2Q3.locomotion)
  output
}

minimum.frames<-function(x) min(apply(X = x,1, length)) #Determine the final frame as there is variation
frame.cut<-function(mat) mat[,1:finalframe]
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

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
Q1Q4.l<-vector("list", length = length(ugenotypes));names(Q1Q4.l)<-ugenotypes
Q2Q3.l<-vector("list", length = length(ugenotypes));names(Q2Q3.l)<-ugenotypes
for(i in 1:length(ugenotypes)) {
  files<-grep(paste0("_", ugenotypes[i]),
              dir, value=TRUE, fixed=TRUE)
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    PI[[i]]<-rbind(PI[[i]], PIextract(data))
    Q1Q4.l[[i]]<-rbind(Q1Q4.l[[i]],Q1Q4.locomotion(data))
    Q2Q3.l[[i]]<-rbind(Q2Q3.l[[i]],Q2Q3.locomotion(data))
  }
}
#Fix the length of the frames for further analysis, just with locomotion for now

finalframe<-min(melt(lapply(Q1Q4.l, FUN = minimum.frames))[,1])
Q1Q4.l<-lapply(Q1Q4.l, frame.cut)

finalframe<-min(melt(lapply(Q2Q3.l, FUN = minimum.frames))[,1])
Q2Q3.l<-lapply(Q2Q3.l, frame.cut)

index<-c(1:finalframe)/30 #This is the framerate from the camerasettings.json file

#Locomoation is calculated as (L-L0)/L0. Note that Q3Q4 stimulation happens first!
single_loco<-function(x) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(x)))
  First30<-colMeans(m[m$seconds<=30 ,])[-1]
  Stim1<-colMeans(m[m$seconds>=30 & m$seconds<=60 ,])[-1]
  Stim2<-colMeans(m[m$seconds>=30 & m$seconds<=60 ,])[-1]
  data.frame(First30, Stim1, Stim2)
}
mat_single_loco<-function(x) apply(x, 1, single_loco) #These function calculate average locomotion proxy for each genotype at stim and control timepoints
Q1Q4single<-melt(sapply(Q1Q4.l, mat_single_loco)) #All the averages for Q1Q4 locomotion at control and stim timepoints
Q1Q4single$quad<-"Q1Q4"
Q2Q3single<-melt(sapply(Q2Q3.l, mat_single_loco)) #All the averages for Q2Q3 locomotion at control and stim timepoints
Q2Q3single$quad<-"Q2Q3"

locomotion<-rbind(Q1Q4single, Q2Q3single)
names(locomotion)<-c("Stim", "value", "N", "Genotype", "quad")
locomotion_singleVal<-data.frame()
for(i in 1:length(ugenotypes)) {
  data<-filter(locomotion, Genotype==ugenotypes[i]) #First select everything with the same genotype

  #calculate L0, all quads and the First30
  L0<-filter(data, Stim=="First30")
  L0<-aggregate(x = L0$value, by=list( L0$N), FUN=mean)
  names(L0)<-c("N", "value")

  #caluclate L, which is the average of the stimulation quadrants' averages
  L<-rbind(filter(data, Stim=="Stim1" & quad=="Q2Q3"), filter(data,Stim=="Stim2" & quad=="Q1Q4"))
  L<-aggregate(x = L$value, by=list( L$N), FUN=mean)
  names(L)<-c("N", "value")

  #Perform the delta L calculation
  deltaL<-((L-L0)/L0)[,-1] #This is a df of all the single value differences in locomotion proxy. One for each experiment.
  deltaL<-data.frame(Genotype=ugenotypes[i], deltaL)
  locomotion_singleVal<-rbind(locomotion_singleVal, deltaL)
}
locomotion_singleVal<-arrange(locomotion_singleVal, desc(Genotype=="EmptySp"))
locomotion_singleVal[,1]<-as.factor(as.character(locomotion_singleVal[,1])) #Reorganise the factor levels as dunn test function takes first level
save(locomotion_singleVal, file=paste0(getwd(),"/locomotion_data.rda"))


#Determine significance and plot the data
#Function to pull out the p-values
pvals<-as.data.frame(dunn.test.control(x = locomotion_singleVal$deltaL,
                                       g= as.factor(locomotion_singleVal$Genotype), p.adjust.method = "fdr")$p.value)
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  locomotion_singleVal, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.1, "Significant", "Not Significant")

g<-ggplot(data=pvals, aes(x=reorder(Genotype, deltaL), y=deltaL))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="deltaL",title="") #Titles
g<-g+theme(legend.title=element_blank())
g


#After running the analyses for all the screen data, let's put everything together for locomotion
#ANALYSIS III: Combining data from both screening sessions, plotting both and analysing repeats
setwd("/Volumes/Data/BehaviourData/")
Dec2015<-loadRData("./Mike_newrig_Dec2015_screen/locomotion_data.rda")
Dec2015<-filter(Dec2015, Genotype!="Empty") #Bug in the DunnTest code, need to remove this to use EmptySp as control
Sept2016<-loadRData("./Mike_newrig_Sept2016_screen/locomotion_data.rda")
all.data<-rbind(Dec2015,Sept2016)
#Remove unwanted samples
all.data<-filter(all.data, Genotype!="test"| Genotype!="MB83c")
all.data<-arrange(all.data, desc(Genotype=="EmptySp"))
all.data[,1]<-as.factor(as.character(all.data[,1])) #Reorganise the factor levels as dunn test function takes first level
pvals<-as.data.frame(dunn.test.control(x = all.data$deltaL,
                                       g= as.factor(all.data$Genotype), p.adjust.method = "fdr")$p.value)
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  all.data, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.1, "Significant", "Not Significant")

g<-ggplot(data=pvals, aes(x=reorder(Genotype, deltaL, median), y=deltaL))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="deltaL",title="") #Titles
g<-g+theme(legend.title=element_blank())
g
ggsave(filename = "allscreen_locomotionproxy_FDR10.pdf", width=16
       ,height =9 ,plot=g, path=".")

