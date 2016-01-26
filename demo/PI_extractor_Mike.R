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

#Plot everything with the empty control, cut the Y-axis at .5
for(i in 2:(length(ugenotypes))) {
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


#Function to calculate the single-value mean from a desired window of a vector. 
#The stimulation windows are 30-60sec and 90-120sec. Yoshi uses last 5 seconds.
singlePI<-function(PIseries, w1s=55, w1e=60, w2s=115, w2e=120) {
  m<-data.frame(cbind(seconds= c(1:3588)/30,as.data.frame(PIseries)))
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
Empty<-filter(all.singlePI.df, Genotype=="Empty")
L989<-filter(all.singlePI.df, Genotype=="L989")
g<-ggplot(data=all.singlePI.df, aes(x=reorder(Genotype, meanPI)
                              , y=meanPI))
g<-g+geom_point(stat = "identity", size=2, alpha=.75)
g<-g+geom_point(data=Empty, color="magenta")
g<-g+geom_point(data=L989, color="green")
g<-g+geom_hline(yintercept = Empty$meanPI)
g<-g+geom_errorbar(aes(ymin=meanPI-sem, ymax=meanPI+sem))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1)) #horizontal text
g






g<-g+geom_point(data=)



