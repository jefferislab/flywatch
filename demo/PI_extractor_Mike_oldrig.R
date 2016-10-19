#This version of the PI extraction  code is optimised for the older rig and in this case is looking 
# specifcally at the iPN line versus empty control 

#Note you will need to exclude data with n=1 and change finalframe variable and index variable 
#if you end up using data from the old rig (different frame rate).
#Set the working directory and load up the required packages and functions
library(tiff)
library(rmngb) 
library(ggplot2)
library(reshape2)
library(dplyr)
library(PMCMR)
library(car)
library(FSA)

PIext<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2))
  PItime<-sapply(data, "[", s)
  PItime
}
#Sort through the directory and pull out the TIFFs.
dir<-list.files(pattern=".tif$", recursive=TRUE)

#Find the different genotypes by pulling out all the genotype data from TIFF filenames
codes<-strsplit(dir, "_")
genotypes<-sapply(codes, "[", 2)
dir2<-data.frame(Tiffs=dir, line=genotypes)
ugenotypes<-unique(genotypes)
print("These are the genotypes to be analysed. Please check for naming errors:");print(ugenotypes)

#Initialise and run the PIext function accross all genotypes. Will produce a list, each containing
#vectors of the PI for each genotype
totalPI<-vector("list", length = length(ugenotypes))
names(totalPI)<-ugenotypes
for(i in 1:length(ugenotypes)) {
  files<-grep(paste0("_", ugenotypes[i]), 
              dir, value=TRUE, fixed=TRUE)
  if (length(files)==1) next
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE) 
    totalPI[[i]]<-rbind(totalPI[[i]], PIext(data))
  }
}

#A function to equalise frames at the end so I can 
# easily make dataframes. The number depends on the arena.  
colGet<-function(mat, finalframe) mat[,1:finalframe]
totalPI<-lapply(totalPI, colGet, finalframe=3642)



#Read out the means of each genotype and give a list 
meanPI<-vector("list", length=length(ugenotypes))
names(meanPI)<-ugenotypes
meanPI<-lapply(totalPI, colMeans)

#Calculate the SEM
varsPI<-lapply(totalPI, colVars)
n<-sapply(sapply(totalPI, subset, select=1),  length)
semPI<-vector("list", length=length(ugenotypes))
names(semPI)<-paste0(ugenotypes, "sem")
for(i in 1:length(ugenotypes)) {
  semPI[[i]]<- sqrt(varsPI[[i]]/n[i])
}

#Convert everything to one df for ggplot2.
index<-c(1:3642)/30 #This is the framerate from the camerasettings.json file
PI_df<-cbind(cbind(seconds= index,as.data.frame(meanPI))
             ,as.data.frame(semPI))

#Plot everything with the empty control, cut the Y-axis at .5
for(i in 2:(length(ugenotypes)+1)) {
  g<-ggplot(data = PI_df,aes(x=seconds))
  g<-g+geom_line(aes(y=PI_df[,i], color=names(PI_df)[i]))
  g<-g+geom_line(aes(y=MB83C, color="MB83C"))
  g<-g+geom_ribbon(aes(ymin=PI_df[,i]-PI_df[,i+length(ugenotypes)], 
                       ymax=PI_df[,i]+PI_df[,i+length(ugenotypes)] )
                   , alpha=.3)
  g<-g+geom_ribbon(aes(ymin=MB83C-MB83Csem, ymax=MB83C+MB83Csem), alpha=.3)
  g<-g+geom_rect(xmax=60, xmin=30, ymax=0, ymin=-1, alpha=0.002, fill="red")
  g<-g+geom_rect(xmax=120, xmin=90, ymax=1, ymin=0, alpha=0.002, fill="red")
  g<-g+coord_cartesian(xlim = c(0, 120), ylim=c(-1,1),expand=FALSE) 
  g<-g+labs(x="Time (Seconds)", y="PI",title="20XUAS-ChrimsonR") 
  g<-g+theme(legend.title=element_blank()) #Turn off the legend title
  filename<-paste0(names(PI_df[i]), "_meanPI.png")
  ggsave(filename = filename, plot=g, path=".")
}


#Function to calculate the single-value mean from a desired window of a vector. 
#The stimulation windows are 30-60sec and 90-120sec. Yoshi uses last 5 seconds.
singlePI<-function(PIseries, w1s=55, w1e=60, w2s=115, w2e=120) {
  m<-data.frame(cbind(seconds= c(1:3642)/30,as.data.frame(PIseries)))
  m1<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
  m2<-colMeans(m[m$seconds>=w2s & m$seconds<=w2e ,])[-1]
  colMeans(rbind(-1*m1,m2))
}

#Function to run singlePI through each row in a matrix. Note first equalise frames using 
#code above. 
mat.singlePI<-function(x) apply(x, 1, singlePI)

#Run mat.singlePI() through all elements of the list 
all.singlePI<-lapply(totalPI, mat.singlePI) 

#Calculate the mean of each genotype under the specified window and get a summary
all.singlePI.means<-sapply(all.singlePI, mean)
summary(all.singlePI.means)

#Calculate the SEM of each genotype under the specified window
sem<-function(x) sqrt(var(x)/length(x))
all.singlePI.sem<-sapply(all.singlePI, sem)

#Convert these to simple dataframes and start plotting. Note to use reorder 
#function to organise the graph as ggplot2 doesn't care about arranged dfs
all.singlePI.df<-data.frame(Genotype=names(all.singlePI.means)
                            ,meanPI=unname(all.singlePI.means)
                            , sem=unname(all.singlePI.sem))

#Remove genotypes you do not want in final analysis
all.singlePI.df_clean<-all.singlePI.df

#Plot cleaned data
g<-ggplot(data=all.singlePI.df_clean, aes(x=reorder(Genotype, meanPI)
                                          , y=meanPI))
g<-g+geom_point(stat = "identity", size=3)
g<-g+geom_hline(yintercept = Empty$meanPI)
g<-g+geom_errorbar(aes(ymin=meanPI-sem, ymax=meanPI+sem))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1)) #horizontal text
g<-g+labs(x="Genotype", y="mean PI (5sec window) with SEM",title="") #Titles
g

#Turn all.singlePI into a tidy dataframe 
all.singlePI.melt<-melt(all.singlePI[2:4])
names(all.singlePI.melt)<-c("PI", "Genotype")
g<-ggplot(data=all.singlePI.melt, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Genotype))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15)) #horizontal text
g<-g+labs(x="", y="PI (5sec window)")
g<-g+theme(legend.position="none")
g
#Remove unwanted genotypes 
all.singlePI.melt2 <-filter(all.singlePI.melt, Genotype=="MB83C" | Genotype=="SS01113")
g<-ggplot(data=all.singlePI.melt2, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Genotype))
g<-g+theme(axis.ticks = element_blank(), axis.text.x = element_blank())
g<-g+labs(x="", y="Performance Index")
g<-g+coord_cartesian(xlim=c(.5,2.5), ylim=c(-1,1), expand=FALSE) 
g

#Do a wilcoxon test 
all.singlePI.melt_clean<-filter(all.singlePI.melt, Genotype=="MB83C" | Genotype=="SS01113")
all.singlePI.melt_clean$Genotype<-as.factor(all.singlePI.melt_clean$Genotype)
wilcox.test(PI ~ Genotype,data=all.singlePI.melt_clean)

#For the 37G11 data, run a KW test followed by a dunn's test
leveneTest(PI~Genotype, data = all.singlePI.melt) #Test for heteroskedasticity which would violate assumptions
kruskal.test(all.singlePI.melt)  #First a Kruskal-Wallis omnibus test, implying significant differences
dunnTest(PI~Genotype, data = all.singlePI.melt,method="bonferroni") #All-v-All from the FSA package



