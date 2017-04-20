#Code to extract tracking data from the tif output of Yoshi's marco

#Load up required packages and functions
library(tiff)
library(dplyr)
library(reshape2)
library(PMCMR)
library(ggplot2)
library(car)


#Set the control genotype
control<-"EmptySp"
#Extract the different metrics from the tifs
PIextract<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2))
  PItime<-sapply(data, "[", s)
  PItime
}
metric_extract2<- function(data, fly, metric=c("DistCenter", "DistBorder","DistNeigh", "CforLoco", "CTurning")) {
  if(metric=="DistCenter") col<-13
  if(metric=="DistBorder") col<-14
  if(metric=="DistNeigh") col<-26
  if(metric=="CforLoco") col<-30
  if(metric=="CTurning") col<-31
  metmatrix<-matrix(ncol=2, byrow=TRUE, data=c(col,fly))
  quadmatrix<-matrix(ncol=2, byrow=TRUE, data=c(12,fly))
  quadrant<-sapply(data, "[", quadmatrix) #pulls the quadrant position of the fly
  metric<-sapply(data, "[", metmatrix)
  sec<-(1:length(metric))/30
  data.frame(quadrant, metric, sec)
}

#Functions to classify and average the metric across distinct quadrant positions for each fly.
#Used on Mean_Of_Exp.
First30_perfly_mean<-function(df)  {
  df<-filter(df, sec<=30) #Metric for all four quadrants before stimulation
  mean(df$metric, na.rm=TRUE)
} #Mean metric for all quadrants before stimulation
Stim1_perfly_mean<-function(df)  {
  df<-filter(df, sec>30 & sec<=60)
  df<-filter(df, quadrant==2 | quadrant==3)
  mean(df$metric, na.rm=TRUE)
} #Mean Metric for the first pair of stimulation quadrants
Stim2_perfly_mean<-function(df)  {
  df<-filter(df, sec>90)
  df<-filter(df, quadrant==1 | quadrant==4)
  mean(df$metric, na.rm=TRUE)
} #Mean Metric for the first pair of stimulation quadrants
Mean_Of_Exp<-function(metric.list) {
  #Pull the means of the metric from the metric.list for each stimulation period
  #and average by fly and then accross flies
  First30meanperfly<-melt(lapply(X = metric.list, FUN = First30_perfly_mean))
  First30meanperexp<-mean(First30meanperfly$value, na.rm = TRUE)

  Stim1meanperfly<-melt(lapply(X = metric.list, FUN = Stim1_perfly_mean))
  Stim1meanperexp<-mean(Stim1meanperfly$value, na.rm = TRUE)

  Stim2meanperfly<-melt(lapply(X = metric.list, FUN = Stim2_perfly_mean))
  Stim2meanperexp<-mean(Stim2meanperfly$value, na.rm = TRUE)

  data.frame(First30=First30meanperexp, Stim1=Stim1meanperexp, Stim2=Stim2meanperexp)
}
delta_metric<-function(df) {
  df<-mutate(df, M=((Stim1+Stim2)/2))
  df<-mutate(df, deltaM_M=(M-First30)/abs(First30))
  metricSingleVal<-select(df, Genotype, deltaM_M)
  metricSingleVal
} #Calculate the delta metric value.
calculate_significants<-function(dataframe, type=c("baseline", "deltametric"), p=.10){
  #Calculate the significant cell-types and label them in the dataframe
  #Uses a posthoc dunn's control test and 10% FDR
  if(type=="baseline") {
    pvals<-as.data.frame(dunn.test.control(x = dataframe$First30,
                                           g= as.factor(dataframe$Genotype),
                                           p.adjust.method = "fdr")$p.value)
  }
  if(type=="deltametric") {
    pvals<-as.data.frame(dunn.test.control(x = dataframe$deltaM_M,
                                           g= as.factor(dataframe$Genotype),
                                           p.adjust.method = "fdr")$p.value)
  }

  pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
  rownames(pvals)<-NULL #unname function doesn't work here.
  names(pvals)<-c("Genotype", "pvalue_v_Empty")
  pvals<-merge(x = pvals, y =  dataframe, by = "Genotype", all=TRUE)
  pvals$Valence<-ifelse(pvals$pvalue_v_Empty<p, "Significant", "Not Significant")
  pvals
}
signif_boxplot<-function(data=df, type=c("baseline", "deltametric")) {
  #Plot the data as boxplots with the statistical significance

  if(type=="baseline") {
    g<-ggplot(data=data, aes(x=reorder(Genotype, First30), y=First30))
    name<-paste0(metric, "_baseline.pdf")
  }
  if(type=="deltametric") {
    g<-ggplot(data=data, aes(x=reorder(Genotype, deltaM_M), y=deltaM_M))
    name<-paste0(metric, "_delta.pdf")
  }
  g<-g+geom_boxplot(aes(fill=Valence))
  g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
  g<-g+labs(x="Genotype", y="Performance Index",title=metric) #Titles
  g<-g+theme(legend.title=element_blank())
  g
  ggsave(filename = name, width=16
         ,height =9 ,plot=g, path=".")
}
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#First loop through the experiments by genotype
dir<-list.files(pattern="TrackingResults.tif$", recursive=TRUE) #Pull out the correct tifs
genotypes<-sapply(strsplit(dir, "_"), "[", 6)
genotypes<-sapply(strsplit(genotypes, "/"), "[", 1)
ugenotypes<-unique(genotypes)
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes)

#Extract the "per experiment" mean of each metric for all flies in desired quadrant and times.
DistCenter.mean.per.exp<-data.frame()
DistBorder.mean.per.exp<-data.frame()
DistNeigh.mean.per.exp<-data.frame()
CforLoco.mean.per.exp<-data.frame()
CTurning.mean.per.exp<-data.frame()
totalPI<-vector("list", length = length(ugenotypes))
for(i in 1:length(ugenotypes)) {
  #First loop through the genotypes
  files<-grep(paste0("_", ugenotypes[i])
              ,dir, value=TRUE, fixed=TRUE)
  #Loop through each experiment
  for(j in 1:length(files))  {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    MaxN<-max(sapply(data, "[", matrix(ncol=2, byrow=TRUE, data=c(1,1))))
    metric.list.DistCenter<-list()
    metric.list.DistBorder<-list()
    metric.list.DistNeigh<-list()
    metric.list.CforLoco<-list()
    metric.list.CTurning<-list()

    totalPI[[i]]<-rbind(totalPI[[i]], PIextract(data))

    #Loop through each fly to extract the metric we want
    for(k in 1:MaxN)   {
      extract<-metric_extract2(data, k, metric = "DistCenter")
      if(mean(extract$metric)==0) next #remove errors
      metric.list.DistCenter[[k]]<-extract
      names(metric.list.DistCenter)[k]<-paste0("Fly", k)

      extract<-metric_extract2(data, k, metric = "DistBorder")
      if(mean(extract$metric)==0) next #remove errors
      metric.list.DistBorder[[k]]<-extract
      names(metric.list.DistBorder)[k]<-paste0("Fly", k)

      extract<-metric_extract2(data, k, metric = "DistNeigh")
      if(mean(extract$metric)==0) next #remove errors
      metric.list.DistNeigh[[k]]<-extract
      names(metric.list.DistNeigh)[k]<-paste0("Fly", k)

      extract<-metric_extract2(data, k, metric = "CforLoco")
      if(mean(extract$metric)==0) next #remove errors
      metric.list.CforLoco[[k]]<-extract
      names(metric.list.CforLoco)[k]<-paste0("Fly", k)

      extract<-metric_extract2(data, k, metric = "CTurning")
      if(mean(extract$metric)==0) next #remove errors
      metric.list.CTurning[[k]]<-extract
      names(metric.list.CTurning)[k]<-paste0("Fly", k)
    }
    #So now we have metric.list, a list object that has our metric against seconds and quadrants
    #for the whole experiment for all the flies.Rbind the mean to the main output df for this experiment
    outputDistCenter<-Mean_Of_Exp(metric.list.DistCenter)
    outputDistCenter$Genotype<-ugenotypes[i]
    DistCenter.mean.per.exp<-rbind(DistCenter.mean.per.exp, outputDistCenter)

    outputDistBorder<-Mean_Of_Exp(metric.list.DistBorder)
    outputDistBorder$Genotype<-ugenotypes[i]
    DistBorder.mean.per.exp<-rbind(DistBorder.mean.per.exp, outputDistBorder)

    outputDistNeigh<-Mean_Of_Exp(metric.list.DistNeigh)
    outputDistNeigh$Genotype<-ugenotypes[i]
    DistNeigh.mean.per.exp<-rbind(DistNeigh.mean.per.exp, outputDistNeigh)

    outputCforLoco<-Mean_Of_Exp(metric.list.CforLoco)
    outputCforLoco$Genotype<-ugenotypes[i]
    CforLoco.mean.per.exp<-rbind(CforLoco.mean.per.exp, outputCforLoco)

    outputCTurning<-Mean_Of_Exp(metric.list.CTurning)
    outputCTurning$Genotype<-ugenotypes[i]
    CTurning.mean.per.exp<-rbind(CTurning.mean.per.exp, outputCTurning)
  }
}

#PART I: Analyse PI
minimum.frames<-function(x) min(apply(X = x,1, length))
finalframe<-min(melt(lapply(totalPI, FUN = minimum.frames))[,1])
frame.cut<-function(mat) mat[,1:finalframe] #Equalise frames at the end
totalPI<-lapply(totalPI, frame.cut)
index<-c(1:finalframe)/30 #This is the framerate from the camerasettings.json file

singlePI<-function(PIseries, w1s=55, w1e=60, w2s=115, w2e=120) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(PIseries)))
  m1<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
  m2<-colMeans(m[m$seconds>=w2s & m$seconds<=w2e ,])[-1]
  colMeans(rbind(-1*m1,m2))
} #Set the PI analysis window
mat.singlePI<-function(x) apply(x, 1, singlePI) #Function to run singlePI through each row in a matrix. Frames must be equalised.
all.singlePI<-lapply(totalPI, mat.singlePI) #Run mat.singlePI() through all elements of the list
all.singlePI.melt<-melt(all.singlePI)
names(all.singlePI.melt)<-c("PI", "Genotype")
all.singlePI.melt<-arrange(all.singlePI.melt, desc(Genotype=="EmptySp"))
all.singlePI.melt<-filter(all.singlePI.melt,Genotype!="test")
save(all.singlePI.melt, file=paste0(getwd(),"/PI_tracking_macro.rda")

#PART II: Analyse the metric data
#Combine all the mean metrics data into one list
all.metrics<-list(DistCenter.mean.per.exp
                  ,DistBorder.mean.per.exp,DistNeigh.mean.per.exp
                  ,CforLoco.mean.per.exp, CTurning.mean.per.exp)
names(all.metrics)<-c("DistCenter", "DistBorder","DistNeigh", "CforLoco", "CTurning")
save(all.metrics, file = "Tracking_all.metrics.rda")

#For each metric in the all.metrics list, calculate the baseline and delta metric
#,run statistics and plot the graphs.
for(i in 1:length(all.metrics)) {
  metric<-names(all.metrics[i])
  baseline<-select(all.metrics[[i]], Genotype, First30)
  baseline<-arrange(baseline, desc(Genotype=="EmptySp"))
  if(metric=="CTurning") {
    baseline$First30<-abs(baseline$First30)
  }
  #Calculate the significant differences between the baselines and plot
  pvals<-calculate_significants(baseline, type="baseline")
  signif_boxplot(data=pvals,type="baseline" )

  #Calculate the single value deltaM/M of the metric
  deltam<-delta_metric(all.metrics[[i]])
  deltam<-arrange(deltam, desc(Genotype=="EmptySp"))
  if(metric=="CTurning") {
    deltam$deltaM_M<-abs(deltam$deltaM_M)
  }
  #Calculate the significant differences for deltaM and plot
  pvals<-calculate_significants(deltam, type="deltametric")
  signif_boxplot(data=pvals,type="deltametric" )
}

#Combine the two screens for full analysis
setwd("/Volumes/Data/BehaviourData")
Dec2015<-loadRData("/Volumes/Data/BehaviourData/Mike_newrig_Dec2015_screen/Tracking_all.metrics.rda")
Sept2016<-loadRData("/Volumes/Data/BehaviourData/Mike_newrig_Sept2016_screen/Tracking_all.metrics.rda")
merge_metric_lists<-function(x, y) {
  output<-vector("list", 5)
  names(output)<-names(x)
  for(i in 1:length(x)){
    output[[i]]<-rbind(x[[i]], y[[i]])
  }
  output
} #Code to concancenate the two lists
total.metrics<-merge_metric_lists(Dec2015,Sept2016)

#Quality Control of the data
remove_unwanted_lines<-function(df){
 df<-filter(df, Genotype!= "11E08")
 df<-filter(df, Genotype!= "53F04")
 df<-filter(df, Genotype!= "empsp")
 df<-filter(df, Genotype!= "Empty")
 df<-filter(df, Genotype!= "MB83C")
 df<-filter(df, Genotype!= "MB083C")
 df
  }
total.metrics<-lapply(total.metrics, FUN = remove_unwanted_lines)
#Calulate and plot the PI data


#Calculate baseline, delta, stats and plot final graphs
for(i in 1:length(total.metrics)) {
  metric<-names(total.metrics[i])
  baseline<-select(total.metrics[[i]], Genotype, First30)
  baseline<-arrange(baseline, desc(Genotype=="EmptySp"))
  if(metric=="CTurning") {
    baseline$First30<-abs(baseline$First30)
  }
  #Calculate the significant differences between the baselines and plot
  pvals<-calculate_significants(baseline, type="baseline")
  signif_boxplot(data=pvals,type="baseline" )

  #Calculate the single value deltaM/M of the metric
  deltam<-delta_metric(total.metrics[[i]])
  deltam<-arrange(deltam, desc(Genotype=="EmptySp"))
  if(metric=="CTurning") {
    deltam$deltaM_M<-abs(deltam$deltaM_M)
  }
  #Calculate the significant differences for deltaM and plot
  pvals<-calculate_significants(deltam, type="deltametric")
  signif_boxplot(data=pvals,type="deltametric" )
}
save(total.metrics, file = "Tracking_all.metrics.rda")


#Fix the code to pull out PI also, code in significant cell-types

