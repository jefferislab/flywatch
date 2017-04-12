#This is code to analyse output from Yoshi's Tracking macro. This is more complicated because it
#tracks each individual fly and therefore requires a different kind of analysis (although you can
# still pull out the PI with my old code)

#Load up required packages and functions
library(tiff)
library(dplyr)
library(reshape2)
library(PMCMR)
library(ggplot2)
library(car)
#Set the control genotype
control<-"empsp"
#Function to extract a metric for a given fly throughout the exp plus what quadrant it is in.
metric_extract<- function(data, fly) {
  #Here we will need to cycle through the columns
  #(31,1) gives us cturning
  metmatrix<-matrix(ncol=2, byrow=TRUE, data=c(31,fly))
  quadmatrix<-matrix(ncol=2, byrow=TRUE, data=c(12,fly))

  quadrant<-sapply(data, "[", quadmatrix) #pulls the quadrant position of the fly
  metric<-sapply(data, "[", metmatrix)
  sec<-(1:length(metric))/30
  data.frame(quadrant, metric, sec)
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
#Used on metric.list output.
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
#This function will run through the metric.list() object and extract the mean of the metric
# for the whole experiment.
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
}

#First loop through the experiments by genotype
dir<-list.files(pattern="TrackingResults.tif$", recursive=TRUE) #Pull out the correct tifs
genotypes<-sapply(strsplit(dir, "_"), "[", 7)
genotypes<-sapply(strsplit(genotypes, "/"), "[", 1)
ugenotypes<-unique(genotypes)
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes)

#PART I: Code to pull out a single metric for each fly.
#Load data and run our functions for a given metric accross all flies of all experiments for each genotype
metric.mean.per.exp<-data.frame()
for(i in 1:length(ugenotypes)) {
    #First loop through the genotypes
    files<-grep(paste0("_", ugenotypes[i])
                ,dir, value=TRUE, fixed=TRUE)
    #Loop through each experiment
    for(j in 1:length(files))  {
        data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
        MaxN<-max(sapply(data, "[", matrix(ncol=2, byrow=TRUE, data=c(1,1))))
        metric.list<-list()
        #Loop through each fly to extract the metric we want
        for(k in 1:MaxN)   {
            extract<-metric_extract(data, k)
            if(mean(extract$metric)==0) next #remove errors
            metric.list[[k]]<-extract
            names(metric.list)[k]<-paste0("Fly", k)
        }
        #So now we have metric.list, a list object that has our metric against seconds and quadrants
        #for the whole experiment for all the flies.Rbind the mean to the main output df for this experiment
        output<-Mean_Of_Exp(metric.list)
        output$Genotype<-ugenotypes[i]
        metric.mean.per.exp<-rbind(metric.mean.per.exp, output)
    }
}

#Calculate the single value deltaM/M of the metric. FIX THIS, eg values that increase will give
#negative delta F/F!
head(metric.mean.per.exp)
metric.mean.per.exp<-mutate(metric.mean.per.exp, M=((Stim1+Stim2)/2))
metric.mean.per.exp<-mutate(metric.mean.per.exp, deltaM_M=(M-First30)/First30)
metricSingleVal<-select(metric.mean.per.exp, Genotype, deltaM_M)

#Statistics on this melted dataframe
pvals<-as.data.frame(dunn.test.control(x = metricSingleVal$deltaM_M,
                                       g= as.factor(metricSingleVal$Genotype), p.adjust.method = "fdr")$p.value)
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  metricSingleVal, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.1, "Significant", "Not Significant")

#Plot the data as significant or not-significant
g<-ggplot(data=pvals, aes(x=reorder(Genotype, deltaM_M), y=deltaM_M))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="deltaL",title="") #Titles
g<-g+theme(legend.title=element_blank())
g

#PART II: Trialling out a way of pulling out all the metrics together and final output as a list of
#metric.mean.per.exp dataframes.
#Will need to loop through each of the major metrics
DistCenter.mean.per.exp<-data.frame()
DistBorder.mean.per.exp<-data.frame()
DistNeigh.mean.per.exp<-data.frame()
CforLoco.mean.per.exp<-data.frame()
CTurning.mean.per.exp<-data.frame()
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

#Combine all the mean metrics data into one list
all.metrics<-list(DistCenter.mean.per.exp
                  ,DistBorder.mean.per.exp,DistNeigh.mean.per.exp
                  ,CforLoco.mean.per.exp, CTurning.mean.per.exp)
names(all.metrics)<-c("DistCenter", "DistBorder","DistNeigh", "CforLoco", "CTurning")
save(all.metrics, file = "Tracking_all.metrics.rda")

#Will make a for loop to run through each entry and metric. Want both baseline and delta metric w/ stats
#Want a text file with the summary statistical data
#adjust to have Empty-Split First

metric<-names(all.metrics[1])
baseline<-select(all.metrics[[1]], Genotype, First30)
baseline<-arrange(baseline, desc(Genotype=="EmptySp"))

#Add the summary statistics to a text file
sink_stats<-function(path, type=c("baseline", "deltametric")) {
  sink(path, append = TRUE)
  print("/"); print("/")
  metric
  if(type="baseline") {
      leveneTest(First30~Genotype, data = baseline)
      print("---")
      kruskal.test(baseline)
      print("---")
      dunn.test.control(x = baseline$First30, g= as.factor(baseline$Genotype), p.adjust.method = "fdr")

     }
  if(type="deltametric"){ #Fix this and match the variables
      leveneTest(deltaM_M~Genotype, data = )
      print("---")
      kruskal.test(baseline)
      print("---")
      dunn.test.control(x = baseline$First30, g= as.factor(baseline$Genotype), p.adjust.method = "fdr")
  }
  print("---")
  print("/")
  sink()
}
calculate_significants<-function(dataframe, type=c("baseline", "deltametric"), p=.10){
  if(type="baseline") {
    pvals<-as.data.frame(dunn.test.control(x = dataframe$First30,
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

sink_stats(path="BaseLineStats.txt", type="baseline")
pvals<-calculate_significants(baseline, type="baseline")

signif_boxplot<-function(data)
g<-ggplot(data=pvals, aes(x=reorder(Genotype, First30), y=First30))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title=metric) #Titles
g<-g+theme(legend.title=element_blank())
g

#Come back to this code and notes. Also pull out the PI
#Use lapply or a loop to run all these commands on the dfs within the all.metrics list
#First need to examine average baseline accross genotypes for all metrics
#Run statisitics and output as txt file.
#Plot data with ggplot2 and output as pdf (need to compare baselines)

#Then Calculate the single value deltaM/M of the metric (take absolute value of cuTurning)
#Run statisitics and output as txt file.
#Plot data with ggplot2 and output as pdf
delta_metric<-function(df) {
  df<-mutate(df, M=((Stim1+Stim2)/2))
  df<-mutate(df, deltaM_M=(M-First30)/abs(First30))
  metricSingleVal<-select(df, Genotype, deltaM_M)
  metricSingleVal
}

lapply(all.metrics, delta_metric)

dunn.test.control(x = metricSingleVal$deltaM_M,
                  g= as.factor(metricSingleVal$Genotype), p.adjust.method = "fdr")

