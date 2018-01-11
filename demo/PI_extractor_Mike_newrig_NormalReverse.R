#This code averages accross the normal 4 quadrant protocol and a reversed version.
#This allows for us to control for any quadrant preference that might exist.

library(tiff)
library(rmngb)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xlsx)
library(PMCMR)
library(car)

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
minimum.frames<-function(x) min(apply(X = x,1, length))
frame.cut<-function(mat) mat[,1:finalframe]
singlePI<-function(PIseries, w1s=55, w1e=60, w2s=115, w2e=120) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(PIseries)))
  m1<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
  m2<-colMeans(m[m$seconds>=w2s & m$seconds<=w2e ,])[-1]
  colMeans(rbind(-1*m1,m2))
}
mat.singlePI<-function(x) apply(x, 1, singlePI) #Function to run singlePI through each row in a matrix. Frames must be equalised.

#Sort through the directory and pull out the TIFFs.
dir<-list.files(pattern="measurementsXY.tif$", recursive=TRUE) #Pull out the correct tif

genotypes<-sapply(strsplit(dir, "_"), "[", 5)
ugenotypes<-unique(genotypes)
protocols<-sapply(strsplit(dir, "_"), "[", 4)

#Initialise and run the PIext function accross all genotypes. Will produce a list, each containing
#vectors of the PI for each genotype

#Extract PI for normal protocol
totalPI_normal<-vector("list", length = length(ugenotypes))
names(totalPI_normal)<-ugenotypes
for(i in 1:length(ugenotypes)) {
  files<-grep(paste0("_", ugenotypes[i]),
              dir, value=TRUE, fixed=TRUE)
  files<-grep("_LED-choice_", files, value=TRUE, fixed=TRUE)
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    totalPI_normal[[i]]<-rbind(totalPI_normal[[i]], PIextract(data))
  }
}

#Extract PI for reverse protocol
totalPI_reverse<-vector("list", length = length(ugenotypes))
names(totalPI_reverse)<-ugenotypes
for(i in 1:length(ugenotypes)) {
  files<-grep(paste0("_", ugenotypes[i]),
              dir, value=TRUE, fixed=TRUE)
  files<-grep("_LED-choice-rev_", files, value=TRUE, fixed=TRUE)
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE)
    totalPI_reverse[[i]]<-rbind(totalPI_reverse[[i]], PIextract(data))
  }
}

#Finish processing on the timeseries
totalPI_reverse<-totalPI_reverse[-18]#Need to remove n=1, L728 from totalPI_reverse

finalframe_normal<-min(melt(lapply(totalPI_normal, FUN = minimum.frames))[,1])
finalframe_reverse<-min(melt(lapply(totalPI_reverse, FUN = minimum.frames))[,1])

totalPI_normal<-lapply(totalPI_normal, frame.cut)
totalPI_reverse<-lapply(totalPI_reverse, frame.cut)

finalframe<-min(c(finalframe_normal,finalframe_reverse))
index<-c(1:finalframe)/30 #This is the framerate from the camerasettings.json file

#Extract the PI from both sets of data
all.singlePI_normal<-melt(lapply(totalPI_normal, mat.singlePI))
all.singlePI_reverse<-melt(lapply(totalPI_reverse, mat.singlePI))
all.singlePI_normal$value<-all.singlePI_normal$value*(-1)
all.singlePI.melt<-rbind(all.singlePI_normal, all.singlePI_reverse)
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
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.20, "Significant", "Not Significant")
#Plot the data as boxplots and colour by significance (according to the FDR)
g<-ggplot(data=pvals, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme(legend.title=element_blank())
g
ggsave(filename = "allPI_FDR20.pdf",
       width=16 ,height =9, plot=g, path=".")