#     To do:
# 1. Make sure that the tiffextract function can actually pull out the tifs for both rigs
#in case they are called something different.
# 5. Make it possibe to pull out the traces too, not just singe value PI
# 6. Normalise PI to pre olfactory stimulation, to account for place pref

#Code to extract the PIs from the US substitution experiments, run this code on each rig data separately.
#Note that each rig has 2 arenas, ie. two cameras. The naming system on the rigs has each rig
#labelled as arena1 and arena2 which is confusing!

library(tiff)
library(rmngb)
library(ggplot2)
library(reshape2)
library(dplyr)
library(xlsx)
library(PMCMR)
library(car)

#Load up the required functions
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
tiffextract<- function(path) {
  output<-list.files(path = path, pattern="measurementsXY.tif$", recursive =TRUE)
  output<-grep(pattern = "movie_Train1", x = output
               ,invert=TRUE, value = TRUE, fixed = TRUE)
}
PIextract<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2))
  PItime<-sapply(data, "[", s)
  PItime
}
timestamp<-function(file) {
  output<-sapply(strsplit(file, "/"), "[", 2)
  output<-sapply(strsplit(file, "_"), "[", 1)
  output<-sapply(strsplit(output, "/"), "[", 2)
  output
}
metaextract<-function(file) {
  timestamp<-timestamp(file)
  cam<-sapply(strsplit(file, "_"), "[", 3)
  odors<-sapply(strsplit(file, "_"), "[", 4)
  paste0(timestamp,"_", cam,"_", odors)
} #Pulls out details about the arena, odor and cam.
PI.by.Cam.odor<-function(df) {

  #For Cam0
  Cam0.OCT<-filter(df,  Cam=="Cam0" & Odor== "OCT+")
  Cam0.MCH<-filter(df,  Cam=="Cam0" & Odor== "MCH+")
  Cam0.MCH<-transform(Cam0.MCH, value=-value)
  if(nrow(Cam0.MCH)!=nrow(Cam0.OCT)) stop("Incorrect Number of Reciprocals")
  #Finally combine each reciprocal odor trial
  Cam0<-as.data.frame(cbind((Cam0.MCH$value+Cam0.OCT$value)/2,Cam0.OCT$L1))

  #For Cam1
  Cam1.OCT<-filter(df, Cam=="Cam1" & Odor== "OCT+")
  Cam1.MCH<-filter(df, Cam=="Cam1" & Odor== "MCH+")
  Cam1.MCH<-transform(Cam1.MCH, value=-value)
  if(nrow(Cam1.MCH)!=nrow(Cam1.OCT)) stop("Incorrect Number of Reciprocals")
  Cam1<-as.data.frame(cbind((Cam1.MCH$value+Cam1.OCT$value)/2, Cam1.OCT$L1))

  #Output the n, where each n is a reciprocal
  all.singlePI.USsub<-rbind(Cam0, Cam1)
  names(all.singlePI.USsub)<-c("PI", "Genotype")
  all.singlePI.USsub
}


#Get the genotypes
dir<-tiffextract(path=".") #Pull out the correct tifs
genotypes<-sapply(strsplit(dir, "_"), "[", 5) #May need to change the number depending on rig
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
  rownames(totalPI[[i]])<-metaextract(files)
}

#Code to determine the final frame as there is variablity
#in the number of frames captured by the camera and equalise frames
minimum.frames<-function(x) min(apply(X = x,1, length))
finalframe<-min(melt(lapply(totalPI, FUN = minimum.frames))[,1])
frame.cut<-function(mat) mat[,1:finalframe]
totalPI<-lapply(totalPI, frame.cut)

#ANALYSIS of US-substitution by pulling out the PI from the last 30seconds (same as Yoshi uses)
singlePI<-function(PIseries, w1s=90, w1e=120) {
  m<-data.frame(cbind(seconds= c(1:finalframe)/30,as.data.frame(PIseries)))
  m<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
}
mat.singlePI<-function(x) apply(x, 1, singlePI) #Function to run singlePI through each row in a matrix. Frames must be equalised.
all.singlePI<-lapply(totalPI, mat.singlePI) #Run mat.singlePI() through all elements of the list
all.singlePI<-lapply(all.singlePI, as.matrix) #Make each list element a matrix so rownames are kept by melt
all.singlePI.melt<-melt(all.singlePI)
all.singlePI.melt<-filter(all.singlePI.melt, L1!="test") #Remove test recordings

#Create new df with separate details
Time<-sapply(strsplit(as.character(all.singlePI.melt$Var1), "_"), "[", 1)
Cam<-sapply(strsplit(as.character(all.singlePI.melt$Var1), "_"), "[", 2)
Odor<-sapply(strsplit(as.character(all.singlePI.melt$Var1), "_"), "[", 3)
all.singlePI.melt<-cbind(all.singlePI.melt, Time, Cam, Odor)
all.singlePI.melt<-select(all.singlePI.melt, -Var1, -Var2)
all.singlePI.melt<-arrange(all.singlePI.melt, Time)

#Average the reciprocal experiments
final.data<-PI.by.Cam.odor(all.singlePI.melt)
save(final.data, file = "US_sub.rda")

#Run statistics and plot the data
setwd("/Volumes/Samsung_T3/USsubstitution_Jan2017") #Go to the directory with LHS_Rig and RHS_rig
alldata<-rbind(loadRData(fileName = "./LHS_rig/US_sub.rda"),loadRData(fileName = "./RHS_rig/US_sub.rda"))
alldata$PI<-as.numeric(as.character(alldata$PI))
dunn.test.control(x = alldata$PI, g= as.factor(alldata$Genotype), p.adjust.method = "fdr")
#Plot the data
g<-ggplot(alldata, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot()
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme(legend.title=element_blank())
g<-g+geom_point()
g
ggsave(filename = "USsubstitution_PI_boxplot.pdf", plot=g, path=".")
