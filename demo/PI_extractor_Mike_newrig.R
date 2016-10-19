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
totalPI<-lapply(totalPI, colGet, finalframe=3601)

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
index<-c(1:3601)/30 #This is the framerate from the camerasettings.json file
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

#Plot everything with the empty control, 
for(i in 2:(length(ugenotypes)+1)) {
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
  m<-data.frame(cbind(seconds= c(1:3601)/30,as.data.frame(PIseries)))
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


#To do: remove Paavo's stuff, MB83C errors, empsps. Do this properly with grep
all.singlePI.df_clean<-all.singlePI.df[(-1)*c(2, 5:8, 18, 19, 50, 46,33),]

#Plot cleaned data
g<-ggplot(data=all.singlePI.df_clean, aes(x=reorder(Genotype, meanPI)
                                    , y=meanPI))
g<-g+geom_point(stat = "identity", size=3)
g<-g+geom_hline(yintercept = Empty$meanPI)
g<-g+geom_errorbar(aes(ymin=meanPI-sem, ymax=meanPI+sem))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1)) #horizontal text
g<-g+labs(x="Genotype", y="mean PI (5sec window) with SEM",title="20XUAS-ChrimsonR") #Titles
g

#Turn all.singlePI into a tidy dataframe and plot the data as boxplots
all.singlePI.melt<-melt(all.singlePI)
names(all.singlePI.melt)<-c("PI", "Genotype")
g<-ggplot(data=all.singlePI.melt, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Genotype))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1)) #horizontal text
g<-g+labs(x="Genotype", y="mean PI (5sec window) with SEM",title="20XUAS-ChrimsonR") #Titles
g<-g+theme(legend.position="none")
g

#Remove Paavo's GAL4 lines from my analysis, melt and rename the df
all.singlePI<-all.singlePI[(-1)*c(2, 5:8, 18, 19, 50, 46,33)] 
all.singlePI.melt<-melt(all.singlePI)
names(all.singlePI.melt)<-c("PI", "Genotype")

#First run KW test. Then most people run loads of Mann-Whitney tests and correct the p-values. 
#This is not ideal because the ranks used by K-W is not the same as M-W. Use Dunn's test with a control 
# and fdr, instead as a sort of "non-parametric Dunnett's test".

leveneTest(PI~Genotype, data = all.singlePI.melt) #Test for heteroskedasticity which would violate assumptions
kruskal.test(all.singlePI)  #First a Kruskal-Wallis omnibus test, implying significant differences
dunn.test.control(x = all.singlePI.melt$PI, g= as.factor(all.singlePI.melt$Genotype), p.adjust.method = "fdr") 
pvals<-as.data.frame(dunn.test.control(x = all.singlePI.melt$PI,
                                       g= as.factor(all.singlePI.melt$Genotype), p.adjust.method = "fdr")$p.value)
table(pvals<0.05) #Print out how many statistically significant differences we found 

pvals<-cbind(dimnames(pvals)[[1]],as.data.frame(pvals)) #Some fudging to switch the genotypes to a column
rownames(pvals)<-NULL #unname function doesn't work here. 
names(pvals)<-c("Genotype", "pvalue_v_Empty")
pvals<-merge(x = pvals, y =  all.singlePI.melt, by = "Genotype", all=TRUE)
pvals$Valence<-ifelse(pvals$pvalue_v_Empty<0.05, "Significant", "Not Significant")

#Plot these statistics onto our graph
g<-ggplot(data=pvals, aes(x=reorder(Genotype, PI), y=PI))
g<-g+geom_boxplot(aes(fill=Valence))
g<-g+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) #horizontal text
g<-g+labs(x="Genotype", y="Performance Index",title="") #Titles
g<-g+theme(legend.title=element_blank())
g


#Examine and plot the distribution of the clusters/cell-types examined
LineSum<-read.xlsx(file = "Line_Summary.xlsx", sheetIndex = 1)
LineSum.tab<-as.data.frame(table(LineSum$Splitlines_cluster..Cluster))
names(LineSum.tab)<-c("Cell-Type", "Frequency")
LineSum.tab<-arrange(LineSum.tab, Frequency)
barplot(height=LineSum.tab[,2], names.arg=LineSum.tab[,1], col="steel blue",cex.names=.6,las=2
        , ylab="Frequency")

#Plot the lines and cell-types that were statistically significant 
hits<-melt(unique(as.character(pvals[pvals$pvalue_v_Empty<.05,]$Genotype))[-20])
names(hits)<-"LineCode"
hits<-merge(x = hits, y = LineSum[,c(1, 5)], by="LineCode", all.x=TRUE
            ,all.y = FALSE)
hits<-arrange(.data = hits, Splitlines_cluster..Cluster) #No multiple cell types
hits$signif<-1
hits<-hits[c(-18,-19),] #Remove the controls
names(hits)<-c("LineCode", "Cell-Type", "signif")
hits2<-merge(x = LineSum.tab, y = hits, by = "Cell-Type", all.x = TRUE)
hits2<-arrange(hits2, Frequency)
barplot(height=hits2[,2], names.arg=hits2[,1], col="steel blue",cex.names=.6,las=2
        , ylab="Frequency")
barplot(height=hits2[,4], col="orange",cex.names=.6,las=2
        , ylab="Frequency", add=TRUE)

