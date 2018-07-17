#Code for analysing the Flybowl data with the tracking macro FlyBowlTracking50.ijm
#The framerate for flybowl is ~30
#Load up required packages, locations and functions
library(tiff)
library(dplyr)
library(reshape2)
library(PMCMR)
library(ggplot2)
library(car)
library(matrixStats)
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
metric.extract<- function(data, fly, metric=c("DistCenter", "DistBorder","DistNeigh", "CforLoco", "CTurning")) {
  if(metric=="DistNeigh") col<-26 #This is correct for flybowl
  if(metric=="CforLoco") col<-30
  if(metric=="CTurning") col<-31
  metmatrix<-matrix(ncol=2, byrow=TRUE, data=c(col,fly))
  metric<-sapply(data, "[", metmatrix)
  data.frame( metric)
}
minimum.frames.list<-function(list) {
  #Takes an element of each condition list (which is itself a list of flies) and get min for that condition
  min(unname(unlist(lapply(X=list,FUN = nrow))))
}
frame.cut.list<-function(list) {
  output<-list()
  for(i in 1:length(list)) {
    output[[i]]<-data.frame(metric=list[[i]][1:finalframe,])
  }
  output
}
Cut_final_frame<-function(Total) {
  finalframe<-min(unname(sapply(X = Total, FUN = minimum.frames.list)))
  Total.cut<-lapply(Total, frame.cut.list)
  Total.cut.df<-sapply(X = Total.cut, FUN = as.data.frame)
  Total.cut.df
}
extract_mean_metric<-function(Total.cut.df) {
  mean<-sapply(Total.cut.df, rowMeans)
  mean.melt<-melt(mean)
  time<-c(1:finalframe)/30
  mean.melt<-cbind(time, mean.melt)
  mean.melt<-dplyr::select(mean.melt, time, Var2, value)
  #mean.melt<-filter(mean.melt,value!=0)
  mean.melt
}
RowSEM<-function(df) {
  n<-ncol(df)
  rowVar<-apply(df, 1, var)
  rowSEM<-sqrt(rowVar/n)
  rowSEM
}
extract_SEM_metric<-function(Total.cut.df) {
  sem<-sapply(Total.cut.df, RowSEM)
  sem.melt<-melt(sem)
  time<-c(1:finalframe)/30
  sem.melt<-cbind(time, sem.melt)
  sem.melt<-dplyr::select(sem.melt, time, Var2, value)
  sem.melt
}

#Summarize and organize the experiments by genotype
all_results<-list.files(pattern="TrackingResults.tif$", recursive=TRUE) #Pull out the correct tifs
genotypes<-sapply(strsplit(all_results, "_"), "[", 2)
genotypes<-sapply(strsplit(genotypes, "/"), "[", 2)
ugenotypes<-unique(genotypes)

#Pull out the locomotion metric
FLoco<-list()
for(i in 1:length(ugenotypes)) {
  files_to_analyze<-grep(pattern = ugenotypes[i], x = all_results, value = TRUE)
  temp_Loco<-list()

  #Loop through each experiment
  for(j in 1:length(files_to_analyze))  {
    data<-readTIFF(source = files_to_analyze[j], all=TRUE, as.is=FALSE)
    MaxN<-max(sapply(data, "[", matrix(ncol=2, byrow=TRUE, data=c(1,1)))) #calculate number of flies
    temp_list_Loco<-list()

    #Loop through each fly to extract the metric
    for(k in 1:MaxN)   {
      extract_Loco<-metric.extract(data, k, metric = "CforLoco")
      if(mean(extract_Loco$metric)==0) next #remove errors
      temp_list_Loco[[k]]<-extract_Loco
      names(temp_list_Loco)[k]<-paste0("Fly", k)
    }
    temp_Loco<-c(temp_Loco,temp_list_Loco)
  }
  #Now name this list and concat into a larger list containing the metric list for each condition

  FLoco[[i]]<-temp_Loco
  names(FLoco)[i]<-ugenotypes[i]
}
summary(FLoco)

#Cut the video frames to be the same length and compress each condition to a df
finalframe<-min(sapply(X = FLoco, FUN = minimum.frames.list))
FLoco.cut<-Cut_final_frame(FLoco)

#Plot1: Calculate the mean and SEM for time series plot
mean<-extract_mean_metric(FLoco.cut)
sem<-extract_SEM_metric(FLoco.cut)

#Pick genotypes to plot
mean<-filter(mean, Var2== "empsp" | Var2=="L1129")
sem<-filter(sem, Var2== "empsp" | Var2=="L1129")
names(sem)<-c("time", "Var2", "sem")
sem<-cbind(mean, sem)

g<-ggplot(data = mean, aes(x=time, y=value, colour=Var2))
g<-g+geom_line()
g<-g+geom_ribbon(data=sem, aes(x=time, ymin=value-sem, ymax=value+sem), colour="black", fill="black" , alpha=0.2)
g






