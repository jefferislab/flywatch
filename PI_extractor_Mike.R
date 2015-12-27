#Note you will need to exclude data with n=1 and change finalframe variable and index variable 
#if you end up using data from the old rig (different frame rate).
#Set the working directory and load up the required packages and functions
library(tiff)
library(rmngb) 
library(ggplot2)
library(reshape2)

PIext<- function(data) {
  s<-matrix(ncol=2, byrow=TRUE, data=c(1,2))
  PItime<-sapply(data, "[", s)
  PItime
}
#Sort through the directory and pull out the TIFFs.
dir<-list.files(pattern=".tif$", recursive=TRUE)

#Find the different genotypes by pulling out all the genotype data from TIFF filenames
codes<-strsplit(dir, "_")
genotypes<-sapply(codes, "[", 6)
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
  for(j in 1:length(files)) {
    data<-readTIFF(source = files[j], all=TRUE, as.is=FALSE) 
    totalPI[[i]]<-rbind(totalPI[[i]], PIext(data))
  }
}

#A function to equalise frames at the end so I can 
# easily make dataframes. The number depends on the arena.  
colGet<-function(mat, finalframe) mat[,1:finalframe]
totalPI<-lapply(totalPI, colGet, finalframe=3588)

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

#Plot just the means of all genotypes together on the same plot using ggplot2
index<-c(1:3588)/30 #This is the framerate from the camerasettings.json file
meltPI<-melt(meanPI)
meltPI<-cbind(meltPI, index)
names(meltPI)<-c("PI", "Genotype", "seconds")
g<-ggplot(data = meltPI, aes(x = seconds))
g<-g+geom_line(mapping = aes(y=PI, color=Genotype))
g<-g+geom_rect(xmax=60, xmin=30, ymax=0, ymin=-1, alpha=0.002, fill="red")
g<-g+geom_rect(xmax=120, xmin=90, ymax=1, ymin=0, alpha=0.002, fill="red")
g<-g+coord_cartesian(xlim = c(0, 120), ylim=c(-1,1), expand=FALSE) 
g<-g+labs(x="Time (Seconds)", y="PI",title="20XUAS-ChrimsonR") #Titles
ggsave(filename = "mean_allgenotypes.png", g, path=".")

#Convert everything to one df for ggplot2.
PI_df<-cbind(cbind(seconds= index,as.data.frame(meanPI))
             ,as.data.frame(semPI))

#Plot everything with the empty control.
for(i in 1:length(ugenotypes)+1) {
  g<-ggplot(data = PI_df,aes(x=seconds))
  g<-g+geom_line(aes(y=PI_df[,i], color=names(PI_df)[i]))
  g<-g+geom_line(aes(y=Empty, color="Empty"))
  g<-g+geom_ribbon(aes(ymin=PI_df[,i]-PI_df[,i+length(ugenotypes)], 
                       ymax=PI_df[,i]+PI_df[,i+length(ugenotypes)] )
                   , alpha=.3)
  g<-g+geom_ribbon(aes(ymin=Empty-Emptysem, ymax=Empty+Emptysem), alpha=.3)
  g<-g+geom_rect(xmax=60, xmin=30, ymax=0, ymin=-1, alpha=0.002, fill="red")
  g<-g+geom_rect(xmax=120, xmin=90, ymax=1, ymin=0, alpha=0.002, fill="red")
  g<-g+coord_cartesian(xlim = c(0, 120), ylim=c(-1,1),expand=FALSE) 
  g<-g+labs(x="Time (Seconds)", y="PI",title="20XUAS-ChrimsonR") 
  g<-g+theme(legend.title=element_blank()) #Turn off the legend title
  filename<-paste0(names(PI_df[i]), "_meanPI.png")
  ggsave(filename = filename, plot=g, path=".")
}

#Single PI value, tweak the window width to calculate an appropriate 
# PI. Yoshi uses last 5sec. 
singlePI<-function(df=meanPI, w1s=45, w1e=60, w2s=105, w2e=120) {
        m<-data.frame(cbind(seconds= index,as.data.frame(meanPI)))
        m1<-colMeans(m[m$seconds>=w1s & m$seconds<=w1e ,])[-1]
        m2<-colMeans(m[m$seconds>=w2s & m$seconds<=w2e ,])[-1]
        colMeans(rbind(-1*m1,m2))
}
singlePI()



