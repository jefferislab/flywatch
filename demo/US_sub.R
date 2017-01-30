#Code to extract the PIs from the US substitution experiments
library(tiff)
library(rmngb)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xlsx)
library(PMCMR)
library(car)

#MORE CHALLENGING THAN THOUGHT AS NEED TO COMPARE RIG BY RIG THE DATA
#Load up the required functions
PIextract<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2))
  PItime<-sapply(data, "[", s)
  PItime
}

#Get the genotypes. CHANGE THE NAMES TO GET SECOND VIDEO
dir<-list.files(pattern="\\+_cam_0_measurementsXY.tif$", recursive=TRUE,) #Pull out the correct tif
genotypes<-sapply(strsplit(dir, "_"), "[", 5)
odors<-sapply(strsplit(dir, "_"), "[", 4)
ugenotypes<-unique(genotypes)
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

#ANALYSIS of US-substitution. NOTE FOR US SUB: that we want the second video and the last 60 seconds
#Will need to code in something about the odor used here to properly denote the
singlePI<-function(PIseries, w1s=100, w1e=120) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(PIseries)))
  m1<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
  m1<-(-1*m1) #Depends on the odor used if positive or negative PI
}
mat.singlePI<-function(x) apply(x, 1, singlePI) #Function to run singlePI through each row in a matrix. Frames must be equalised.
all.singlePI<-lapply(totalPI, mat.singlePI) #Run mat.singlePI() through all elements of the list
all.singlePI.melt<-melt(all.singlePI)
names(all.singlePI.melt)<-c("PI", "Genotype")
all.singlePI.melt$PI<-abs(all.singlePI.melt$PI)
save(all.singlePI.melt, file=paste0(getwd(),"/data.rda"))
#Statisics
leveneTest(PI~Genotype, data = all.singlePI.melt) #Test for heteroskedasticity which would violate assumptions
kruskal.test(all.singlePI.melt)  #First a Kruskal-Wallis omnibus test, implying significant differences
dunn.test.control(x = all.singlePI.melt$PI, g= as.factor(all.singlePI.melt$Genotype), p.adjust.method = "fdr")
g<-ggplot(data=all.singlePI.melt, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Genotype))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme(legend.title=element_blank())
g<-g+geom_point(aes(fill=Genotype), pch=21,color="black",  size=3 )
g
ggsave(filename = "PI_boxplot.pdf", plot=g, path=".")

#ANALYSIS II: Plot by odors:
#Get the genotypes
dir<-list.files(pattern="\\+_cam_0_measurementsXY.tif$", recursive=TRUE,) #Pull out the correct tif
genotypes<-sapply(strsplit(dir, "_"), "[", 5)
odors<-sapply(strsplit(dir, "_"), "[", 4)
genotypes.odors<-paste0(genotypes,"_", odors)
ugenotypes.odors<-unique(genotypes.odors)
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes.odors)

#Initialise and run the PIext function accross all genotypes. Will produce a list, each containing
#vectors of the PI for each genotype
totalPI<-vector("list", length = length(genotypes.odors))
names(totalPI)<-genotypes.odors
for(i in 1:length(genotypes.odors)) {
  files<-grep(paste0("_", genotypes.odors[i]),
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







