#This is code to analyse output from Yoshi's Tracking macro. This is more complicated because it
#tracks each individual fly and therefore requires a different kind of analysis (although you can
# still pull out the PI with my old code)

#Load up required packages and functions
library(tiff)
library(dplyr)
library(reshape2)
library(PMCMR)
library(ggplot2)


#Function to extract a metric for a given fly throughout the exp
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

#First loop through the experiments by genotype
dir<-list.files(pattern="TrackingResults.tif$", recursive=TRUE) #Pull out the correct tifs
genotypes<-sapply(strsplit(dir, "_"), "[", 7)
genotypes<-sapply(strsplit(genotypes, "/"), "[", 1)
ugenotypes<-unique(genotypes)
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes)

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
        #for the whole experiment for all the flies.rbind it to the main output df for this experiment
        output<-Mean_Of_Exp(metric.list)
        output$Genotype<-ugenotypes[i]
        metric.mean.per.exp<-rbind(metric.mean.per.exp, output)
    }
}

#Calculate the single value deltaM/M of the metric
head(metric.mean.per.exp)
metric.mean.per.exp<-mutate(metric.mean.per.exp, M=mean(c(Stim1, Stim2)))
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


#Trialling out a way of pulling out all the metrics together and final output as a list of
#metric.mean.per.exp dataframes. Once analysis segment is done, can run all these in parallel
metric_extract<- function(data, fly, metric=c("DistCenter", "DistBorder"
                                              ,"DistNeigh", "CforLoco", "CTurning")) {
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

#Will need to loop through and create a function to automate the plotting.
