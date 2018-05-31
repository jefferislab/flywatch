#Chunk of code to plot behaviour traces more for Figure7 of split-GAL4 paper
g<-ggplot(data = PI_df,aes(x=seconds))
g<-g+geom_line(aes(y=PI_df[,i], color=names(PI_df)[i]))
g<-g+geom_line(aes(y=EmptySp, color="EmptySp"))
g<-g+geom_ribbon(aes(ymin=PI_df[,i]-PI_df[,i+length(ugenotypes)],
                     ymax=PI_df[,i]+PI_df[,i+length(ugenotypes)] )
                 , alpha=.3)
g<-g+geom_ribbon(aes(ymin=EmptySp-EmptySpsem, ymax=EmptySp+EmptySpsem), alpha=.3)
g<-g+geom_rect(xmax=60, xmin=30, ymax=0, ymin=-1,fill="red", alpha=0.002 )
g<-g+geom_rect(xmax=120, xmin=90, ymax=1, ymin=0, , fill="red", alpha=0.002)
g<-g+coord_cartesian(xlim = c(0, 120), ylim=c(-1,1),expand=FALSE)
g<-g+labs(x="Time (Seconds)", y="Performance Index")
g<-g+theme(legend.title=element_blank()) #Turn off the legend title
g<-g+geom_hline(yintercept = 0, alpha=0.5)
g<-g+theme(axis.text.x = element_text(size=12))
g<-g+theme(axis.text.y = element_text(size=12))
g<-g+scale_color_manual(labels = c("Empty Split-GAL4", "L728"), values=c("orange", "steel blue"))
g<-g+theme(legend.text=element_text(size=15))
g+theme(axis.title.x=element_text(size=40))
g+theme(axis.title.y=element_text(size=40))
g