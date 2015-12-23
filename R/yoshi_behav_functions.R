# Greg's little attempt a processing yoshi's behaviour tiffs

#' read xy data from a single behaviour summary file
read_ybr_xy<-function(f) {
  rawd=tiff::readTIFF(f, all=TRUE, as.is=FALSE)
  lapply(rawd, function(x) t(x[2:3,which(x[2,]>0)]))
}

#' Read the summary values (inc PI) from first row of a
#' yoshi behaviour tiff and return a multi timeseries
#' object
#'
#' @examples
#' some_summ=read_ybr_summary("some.tif")
#' plot(some_summ[,c(1,2,7)])
read_ybr_summary<-function(f, ncols=26, start=0, deltat=1/30, ...) {
  rawd=tiff::readTIFF(f, all=TRUE, as.is=FALSE)
  d=as.data.frame(t(sapply(rawd, "[",1, 1:ncols)))
  colnames(d)[1]="nResults"
  colnames(d)[2]="PIn"
  colnames(d)[3:6]=paste0("NQ",1:4)
  colnames(d)[7]="PIpixel"
  colnames(d)[8:11]=paste0("Q",1:4)
  colnames(d)[12:15]=paste0("DQ",1:4)
  colnames(d)[16]="RelativeDQ1DQ4"
  colnames(d)[17]="RelativeDQ2DQ3"
  colnames(d)[18]="DistanceFromCenter.mu"
  colnames(d)[19]="DistanceFromCenter.sd"
  colnames(d)[20]="DistanceFromCenterQ1Q4.mu"
  colnames(d)[21]="DistanceFromCenterQ1Q4.sd"
  colnames(d)[22]="DistanceFromCenterQ2Q3.mu"
  colnames(d)[23]="DistanceFromCenterQ2Q3.sd"
  colnames(d)[24]="DistanceFromBorder.mu"
  colnames(d)[25]="DistanceFromBorder.sd"
  colnames(d)[26]="FractionOfFliesinChoiceZone"
  ts(d, start=start, deltat=deltat, ...)
}

#' calculate the raw displacements between nearest neighbour flies
#' in sequential
#'
calc_raw_displacements<-function(xy){
  distres=list()
  for(f in seq_along(xy)[-1]){
    # find difference in position
    nn=nabor::knn(data = xy[[f]], query=xy[[f-1]], k = 1)
    distres[[f-1]]=c(nn$nn.dists)
  }
  distres
}

median_displacement<-function(xy, deltat=1/30){
  dd=calc_raw_displacements(xy)
  md=sapply(dd, median)
  ts(md,start = 0, deltat = deltat)
}

plot_smoothed_displacement<-function(xy, filterwidth=1, lights=c(30,60,90,120),
                                     lightcol=rgb(1,0,0,alpha=.3), ...){
  if(is.character(xy)) xy=read_ybr_xy(xy)
  mxy=median_displacement(xy)

  f=rep(deltat(mxy)/filterwidth, filterwidth/deltat(mxy))
  sm_mxy=stats::filter(mxy, f)
  plot(sm_mxy, ...)

  rand_ts=ts(sample(mxy), start=start(mxy), deltat = deltat(mxy))
  lines(stats::filter(rand_ts, f), col='red')
  rect(lights[1], par("usr")[3], lights[2],
       par("usr")[4], col=lightcol, border = NA)
  rect(lights[3], par("usr")[3], lights[4],
       par("usr")[4], col=lightcol, border = NA)

}

summarise_tifs<-function(path="."){
  all_tifs=dir(path=path,pattern="tif$", recursive=T)
  timestamps=stringr::str_extract(dirname(all_tifs),"201[0-9]{5}T[0-9]{6}")
  ptimestamps=lubridate::parse_date_time(timestamps,"YmdHMS")
  geno=stringr::str_match(dirname(all_tifs),"_([^_]+)_20X")[,2]
  data.frame(tif=all_tifs,geno=geno, time=ptimestamps, stringsAsFactors = F)
}
