#Note you will need to exclude data with n=1 and change finalframe variable and index variable
#if you end up using data from the old rig (different frame rate).
#Need to clean up older variables, objects and functions!
#Set the working directory and load up the required packages and functions
library(tiff)
library(rmngb)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xlsx)
library(PMCMR)
library(car)

#Remove Cam1 entries
#files<-list.files(path = ".", recursive = TRUE, include.dirs = TRUE)
#Cam1s<-grep(x = files, pattern = "Cam1", fixed = TRUE, value = TRUE)
#file.remove(Cam1s)

#Load up the required functions
PIextract<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2))
  PItime<-sapply(data, "[", s)
  PItime
}
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#Sort through the directory and pull out the TIFFs.
dir<-list.files(pattern="measurementsXY.tif$", recursive=TRUE) #Pull out the correct tif

#Simple code to standardise any naming errors for MB83C
dir.MB83.fix<-gsub(pattern = "*MB83c*", replacement = "MB83C", x = dir)
file.rename(from = dir, to = dir.MB83.fix)

#Find the different genotypes by pulling out all the genotype data from TIFF filenames
if(length(unique(sapply(strsplit(dir, "_"), "[", 7)))==1) {
  genotypes<-sapply(strsplit(dir, "_"), "[", 6)
  ugenotypes<-unique(genotypes)
  ugenotypes<-ugenotypes[(-1)*c(2, 5:8, 18, 19, 33,46,50)]
}
if(length(unique(sapply(strsplit(dir, "_"), "[", 7)))>1) {
  genotypes<-sapply(strsplit(dir, "_"), "[", 7)
  genotypes<-sapply(strsplit(genotypes, "/"), "[", 1)
  ugenotypes<-unique(genotypes)
}
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes)

#Initialise and run the PIext function accross all genotypes. Will produce a list, each containing
#vectors of the PI for each genotype
totalPI<-vector("list", length = length(ugenotypes))
names(totalPI)<-ugenotypes
for(i in 1:length(ugenotypes)) {
  files<-grep(paste0("_", ugenotypes[i]),
              dir, value=TRUE, fixed=TRUE)
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    totalPI[[i]]<-rbind(totalPI[[i]], PIextract(data))
  }
}

#A function to decide the final frame as there is variablity
#in the number of frames captured by the camera.
minimum.frames<-function(x) min(apply(X = x,1, length))
finalframe<-min(melt(lapply(totalPI, FUN = minimum.frames))[,1])

#Equalise frames at the end
frame.cut<-function(mat) mat[,1:finalframe]
totalPI<-lapply(totalPI, frame.cut)
index<-c(1:finalframe)/30 #This is the framerate from the camerasettings.json file

#ANALYSIS I: Plotting the PI as a timeseries
meanPI<-vector("list", length=length(ugenotypes))
names(meanPI)<-ugenotypes
meanPI<-lapply(totalPI, colMeans)
varsPI<-lapply(totalPI, colVars)
n<-sapply(sapply(totalPI, subset, select=1),  length)
semPI<-vector("list", length=length(ugenotypes))
names(semPI)<-paste0(ugenotypes, "sem")
for(i in 1:length(ugenotypes))  semPI[[i]]<- sqrt(varsPI[[i]]/n[i])
PI_df<-cbind(cbind(seconds= index,as.data.frame(meanPI))
             ,as.data.frame(semPI))
#PI_flip<-cbind(select(PI_df, seconds), PI_df[,-1]*(-1)) #For flipping the PI
#PI_df<-PI_flip
#Plot each line with the empty/emptySp control. Set which.
for(i in 2:(length(ugenotypes)+1)) {
  g<-ggplot(data = PI_df,aes(x=seconds))
  g<-g+geom_line(aes(y=PI_df[,i], color=names(PI_df)[i]))
  g<-g+geom_line(aes(y=empsp, color="empsp"))
  g<-g+geom_ribbon(aes(ymin=PI_df[,i]-PI_df[,i+length(ugenotypes)],
                       ymax=PI_df[,i]+PI_df[,i+length(ugenotypes)] )
                   , alpha=.3)
  g<-g+geom_ribbon(aes(ymin=empsp-empspsem, ymax=empsp+empspsem), alpha=.3)
  g<-g+geom_rect(xmax=60, xmin=30, ymax=0, ymin=-1, alpha=0.002, fill="red")
  g<-g+geom_rect(xmax=120, xmin=90, ymax=1, ymin=0, alpha=0.002, fill="red")
  g<-g+coord_cartesian(xlim = c(0, 120), ylim=c(-1,1),expand=FALSE)
  g<-g+labs(x="Time (Seconds)", y="PI",title="20XUAS-csChrimson")
  g<-g+theme(legend.title=element_blank()) #Turn off the legend title
  filename<-paste0(names(PI_df[i]), "_meanPI.png")
  ggsave(filename = filename, plot=g, path=".")
}

#ANALYSIS II: Plotting the single PI values and their statistics
#Define a function to calculate the single-value mean from a desired window of a vector.
#The stimulation windows are 30-60sec and 90-120sec. Yoshi uses last 5 seconds.
singlePI<-function(PIseries, w1s=55, w1e=60, w2s=115, w2e=120) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(PIseries)))
  m1<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
  m2<-colMeans(m[m$seconds>=w2s & m$seconds<=w2e ,])[-1]
  colMeans(rbind(-1*m1,m2))
}
mat.singlePI<-function(x) apply(x, 1, singlePI) #Function to run singlePI through each row in a matrix. Frames must be equalised.
all.singlePI<-lapply(totalPI, mat.singlePI) #Run mat.singlePI() through all elements of the list
all.singlePI.melt<-melt(all.singlePI)
names(all.singlePI.melt)<-c("PI", "Genotype")
all.singlePI.melt<-arrange(all.singlePI.melt, desc(Genotype=="EmptySp"))
all.singlePI.melt<-filter(all.singlePI.melt,Genotype!="test")
save(all.singlePI.melt, file=paste0(getwd(),"/data.rda"))
#Run the statistics
leveneTest(PI~Genotype, data = all.singlePI.melt) #Test for heteroskedasticity which would violate assumptions
SW.test<-function(x) {
   shapiro.test(x)$p.value
}
aggregate(x=all.singlePI.melt$PI, by = list(all.singlePI.melt$Genotype), FUN=SW.test)
#Function to pull out the p-values
kruskal.test(all.singlePI)  #First a Kruskal-Wallis omnibus test, implying significant differences
dunn.test.control(x = all.singlePI.melt$PI, g= as.factor(all.singlePI.melt$Genotype), p.adjust.method = "fdr")
pvals<-as.data.frame(dunn.test.control(x = all.singlePI.melt$PI,
                                       g= as.factor(all.singlePI.melt$Genotype), p.adjust.method = "fdr")$p.value)
table(pvals<0.05) #Print out how many statistically significant differences we found
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  all.singlePI.melt, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.2, "Significant", "Not Significant")
#Plot the data as boxplots and colour by significance (according to the FDR)
g<-ggplot(data=pvals, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme(legend.title=element_blank())
g
ggsave(filename = "AmbientLightON.pdf",
       width=16 ,height =9, plot=g, path=".")

#ANALYSIS III: Combining data from both screening sessions, plotting both and analysing repeats
setwd("/Volumes/Samsung_T3/Behaviour_data")
Dec2015<-loadRData("/Volumes/Data/BehaviourData/Mike_newrig_Dec2015_screen/data.rda")
Dec2015<-filter(Dec2015, Genotype!="Empty") #Bug in the DunnTest code, need to remove this to use EmptySp as control
Sept2016<-loadRData("/Volumes/Data/BehaviourData/Mike_newrig_Sept2016_screen/data.rda")
all.data<-rbind(Dec2015,Sept2016)
#Remove unwanted samples
all.data<-filter(all.data, Genotype!="test"| Genotype!="MB83c")
all.data<-arrange(all.data, desc(Genotype=="EmptySp"))
#Compare the repeated lines in the two screening sessions
repeats<-c("L235", "L728", "L574", "L421", "L141", "L260", "L159", "L123")
repeat.means.Dec2015<- subset(Dec2015, Genotype %in% repeats)
repeat.means.Dec2015<-aggregate(.data = repeat.means.Dec2015, x = repeat.means.Dec2015$PI
          , by=list(repeat.means.Dec2015$Genotype),FUN =  mean)
repeat.means.Sept2016<-subset(Sept2016, Genotype %in% repeats)
repeat.means.Sept2016<-aggregate(.data = repeat.means.Sept2016, x = repeat.means.Sept2016$PI
                                , by=list(repeat.means.Sept2016$Genotype),FUN =  mean)
t.test(x = repeat.means.Sept2016[,2]-repeat.means.Dec2015[,2]) #No significant difference between the repeats in 2015 and 2016
#Compare and plot EmptySp and Empty
controls<-filter(all.data, Genotype=="Empty" | Genotype=="EmptySp")
kruskal.test(controls)  #Compare the two controls, pretty different
g<-ggplot(data = controls, mapping = aes(x=Genotype, y=PI))
g<-g+geom_boxplot(aes(fill=Genotype))
g
ggsave(filename = "PI_comparision_controls.pdf", plot=g, path=".")
#Statistical tests on all the data, using either EmptySplit or Empty as the control
all.data<-arrange(all.data, desc(Genotype=="EmptySp"))
leveneTest(PI~Genotype, data = all.data) #Test for heteroskedasticity which would violate assumptions
kruskal.test(all.data)  #First a Kruskal-Wallis omnibus test, implying significant differences
dunn.test.control(x = all.data$PI, g= as.factor(all.data$Genotype), p.adjust.method = "fdr")
pvals<-as.data.frame(dunn.test.control(x = all.data$PI,
                                       g= as.factor(all.data$Genotype), p.adjust.method = "fdr")$p.value)
p<-0.10
table(pvals<p) #Print out how many statistically significant differences we found
pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here.
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  all.data, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<p, "Significant", "Not Significant")
#Plot the data as boxplots and colour by significance (according to the FDR)
g<-ggplot(data=pvals, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme(legend.title=element_blank())
g
ggsave(filename = "alldata_EmptySplitcontrol_FDR10.pdf", width=16
       ,height =9 ,plot=g, path=".")

#Plot just the significants and colour by cell-type
sig<-filter(pvals, Valence=="Significant" | Genotype=="EmptySp")
sig<-filter(sig, Genotype!="Gr66a")
sig<-filter(sig, Genotype!="L235")
types<-read.xlsx(file = "Line_Summary.xlsx", sheetIndex = 1)
sig_types<-merge(x = sig, by.x="Genotype", y=types, by.y = "LineCode", all.x = TRUE
                 , all.y=FALSE)

g<-ggplot(data=sig_types, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_hline(yintercept =0, alpha=0.6)
g<-g+geom_boxplot(aes(fill=reorder(Clusters..Cluster, PI)), outlier.shape = NA)

g<-g+geom_boxplot(data=filter(sig_types, Genotype=="EmptySp"), col="red", outlier.shape = NA)
g<-g+geom_jitter(alpha=0.2, width=0.1, size=2 )
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme_bw()
g<-g+theme(legend.title=element_blank())
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+theme(axis.text.x = element_text(size=10),  axis.text.y = element_text(size=10))
g
ggsave(filename = "alldata_EmptySplitcontrol_FDR10_significants.pdf",plot=g, path=".")


#ANALYSIS IV: Significance by line cell-type. FIX THIS!
#Examine and plot the distribution of the clusters/cell-types examined
LineSum<-read.xlsx(file = "Line_Summary.xlsx", sheetIndex = 1)[,c(1,3)]
length(unique(LineSum[,2])) #Get the number of cell-types analysed
LineSum.tab<-as.data.frame(table(LineSum$Clusters..Cluster))
names(LineSum.tab)<-c("Cell-Type", "Frequency")
LineSum.tab<-arrange(LineSum.tab, Frequency)
barplot(height=LineSum.tab[,2], names.arg=LineSum.tab[,1], col="red",cex.names=.6,las=2
        , ylab="Frequency")
#Plot the lines and cell-types that were statistically significant
hits<-melt(unique(as.character(pvals[pvals$pvalue_v_Empty<p,]$Genotype)))
names(hits)<-"LineCode"
hits<-merge(x = hits, y = LineSum, by="LineCode", all.x=TRUE
            ,all.y = FALSE)
hits.table<-as.data.frame(table(hits[,2]))
names(hits.table)<-c("Cell-Type", "FreqHit")
hits.merge<-merge(x=hits.table, y=LineSum.tab, by="Cell-Type", all=TRUE)
hits.merge<-arrange(hits.merge, Frequency)
barplot(height=hits.merge[,3], names.arg=LineSum.tab[,1], col="firebrick1",cex.names=.6,las=2
        , ylab="Frequency")
barplot(height=hits.merge[,2], names.arg=LineSum.tab[,1], col="deepskyblue1",cex.names=.6,las=2
        , ylab="Frequency", add=TRUE)


#ANALYSIS V: Comparing two different protocols

#Sort through the directory and pull out the TIFFs, genotypes and protocols
dir<-list.files(pattern="measurementsXY.tif$", recursive=TRUE) #Pull out the correct tif
genotypes<-sapply(strsplit(dir, "_"), "[", 5)
ugenotypes<-unique(genotypes)[-3]
protocols<-unique(sapply(strsplit(dir, "_"), "[", 4))
#Initialise and run the PIext function accross all genotypes. Will produce a list, each containing
#vectors of the PI for each genotype
totalPI<-vector("list", length = length(ugenotypes)*length(protocols))
names(totalPI)<-as.vector(t(outer(protocols, ugenotypes, paste, sep="_")))
for(i in 1:(length(ugenotypes)*length(protocols))) {
  files<-grep(paste0("_", names(totalPI)[i]),
              dir, value=TRUE, fixed=TRUE)
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    totalPI[[i]]<-rbind(totalPI[[i]], PIextract(data))
  }
}
