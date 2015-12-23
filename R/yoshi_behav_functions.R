# Code to process Yoshi's behaviour TIFFs

#' Read in one of Yoshi's behaviour TIFF summary files
#'
#' @description Low level function to read in raw data from a Yoshi behaviour
#'   summary TIFF file. You can use this to customise how the TIFF file is read
#'   in \emph{or} if you want to extract multiple kinds of information from the
#'   same TIFF file using different functions.
#'
#' @param f Path to a TIFF file on disk \bold{or} the data from a TIFF that has
#'   previously been read into R.
#' @param ... Additional arguments passed to \code{\link[tiff]{readTIFF}}
#' @inheritParams tiff::readTIFF
#' @return list of TIFF slices with extra class \code{ybr_raw}
#' @export
#' @family read_ybr
#' @seealso \code{\link[tiff]{readTIFF}}
read_ybr_tiff<-function(f, all=TRUE, as.is=FALSE, ...) {
  if(is.character(f)) {
    l=suppressWarnings(tiff::readTIFF(f, all=all, as.is=as.is, ...))
    class(l)=c('ybr_raw', class(l))
    l
  } else if(inherits(f, 'ybr_raw')) {
    f
  } else stop("Sorry, I don't recognise this\n", f,
              "\nas a path or in memory tiff data!")
}

#' Read xy positions from a single behaviour summary file
#'
#' @param f Path to a TIFF file on disk \bold{or} the data from a TIFF that has
#'   previously been read into R
#' @param ... Additional arguments passed to \code{\link{read_ybr_tiff}} or
#'   eventually to \code{\link[tiff]{readTIFF}}.
#' @family read_ybr
#' @export
read_ybr_xy<-function(f, ...) {
  if(inherits(f, 'ybr_xy')) return(f)
  rawd=read_ybr_tiff(f, ...)
  l=lapply(rawd, function(x) t(x[2:3,which(x[2,]>0)]))
  class(l)=c("ybr_xy")
  l
}

#' Read the summary values from yoshi behaviour TIFF
#'
#' @param start Time of the first observation
#' @param frequency Frequency (i.e. frame rate) for the behavioural observations
#' @param ... Additional arguments passed to \code{\link[stats]{ts}}.
#' @inheritParams read_ybr_xy
#'
#' @return a multi timeseries (\code{\link[stats]{ts}}) object containing 26
#'   columns reporting frame by frame summary statistics.
#'
#' @examples
#' \dontrun{
#' some_summ=read_ybr_summary("some.tif")
#' plot(some_summ[,c(1,2,7)])
#' }
#' @family read_ybr
#' @seealso \code{\link[stats]{ts}}
#' @export
read_ybr_summary<-function(f, start=0, frequency=30, ...) {
  rawd=read_ybr_tiff(f)
  ncols=26
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
  stats::ts(d, start=start, frequency=frequency, ...)
}

#' Raw displacements between nearest neighbour flies in sequential frames
#'
#' Each frame may have a different number of "flies". So all we do is report for
#' each frame \code{n}  the nearest neighbour distance to objects in frame
#' \code{n+1}. For each frame there will therefore be as many distances are
#' there are recognised flies.
#' @param xy List of xy locations per frame produced by
#'   \code{\link{read_ybr_xy}} \emph{or} raw TIFF data or path to TIFF file.
#' @export
#' @family ybr-displacement
ybr_raw_displacements<-function(xy){
  xy=read_ybr_xy(xy)
  distres=list()
  for(f in seq_along(xy)[-1]){
    # find difference in position
    nn=nabor::knn(data = xy[[f]], query=xy[[f-1]], k = 1)
    distres[[f-1]]=c(nn$nn.dists)
  }
  distres
}

#' Median displacement between nearest neighbour flies in sequential frames
#'
#' @inheritParams ybr_raw_displacements
#' @inheritParams read_ybr_summary
#' @return time series
#' @export
#' @family ybr-displacement
ybr_median_displacement<-function(xy, start=0, frequency=30){
  dd=ybr_raw_displacements(xy)
  md=sapply(dd, median)
  ts(md, start = start, frequency = frequency)
}

#' Plot smoothed displacement estimate
#'
#' @inheritParams ybr_raw_displacements
#' @param filter \emph{Either} the width in seconds of a simple smoothing filter
#'   \emph{or} a filter defined according to \code{\link[stats]{filter}}.
#' @param sides Whether the filter is causal i.e. for past values only (the
#'   default) or centered around lag=0. See \code{\link[stats]{filter}} for
#'   details.
#' @param lights A length 2 or more vector defining the lights on/off times for
#'   the experiment (in seconds).
#' @param lightcol The colour to use to plot the lights on epochs
#' @param ... Additional arguments to \code{plot.ts}
#' @family ybr-displacement
#' @export
#' @seealso \code{\link[stats]{filter}}
#' tiffdf=find_ybr_tiffs(system.file("ybr_tiffs", package='flywatch'))
#' plot_smoothed_displacement(tiffdf$tiff[1])
plot_smoothed_displacement<-function(xy, filter=1, sides=1,
                                     lights=c(on1=30, off1=60, on2=90, off2=120),
                                     lightcol=rgb(1,0,0,alpha=.3), ...){
  mxy=ybr_median_displacement(xy)
  # computer filter if required
  f=if(length(filter)==1) {
    rep(deltat(mxy)/filter, filter/deltat(mxy))
  } else filter

  sm_mxy=stats::filter(mxy, f, sides=sides)
  plot(sm_mxy, ylab='median displacement per frame', ...)

  rand_ts=ts(sample(mxy), start=start(mxy), deltat = deltat(mxy))
  lines(stats::filter(rand_ts, f), col='red')
  onidxs=seq.int(from=1, to=length(lights), by=2)
  rect(lights[onidxs], par("usr")[3], lights[onidxs+1],
       par("usr")[4], col=lightcol, border = NA)
}

#' Find all the Yoshi behaviour tiffs in a directory hierarchy
#'
#' @inheritParams base::list.files
#' @return A data.frame containing information about the identified tiffs
#' @export
#' @family read_ybr
#' @examples
#' find_ybr_tiffs(system.file("ybr_tiffs", package='flywatch'))
#'
#' \dontrun{
#' ## Make some plots for all the experiments in a directory
#' tiffdf=find_ybr_tiffs("/path/to/my/flydata")
#' for(i in 1:nrow(tiffdf)){
#'   intiff=tiffdf[i,'tiff']
#'   outpdf=paste0(basename(dirname(intiff)),'.pdf')
#'   pdf(outpdf)
#'   plot_smoothed_displacement(intiff)
#'   dev.off()
#'   cat(".")
#' }
#' }
find_ybr_tiffs<-function(path=".", recursive=TRUE, full.names=TRUE){
  all_tiffs=dir(path=path, pattern="\\.tif$", recursive=recursive,
                full.names=full.names)
  timestamps=stringr::str_extract(dirname(all_tiffs),"201[0-9]{5}T[0-9]{6}")
  ptimestamps=lubridate::parse_date_time(timestamps,"YmdHMS")
  geno=stringr::str_match(dirname(all_tiffs),"_([^_]+)_20X")[,2]
  data.frame(tiff=all_tiffs, geno=geno, time=ptimestamps, stringsAsFactors = F)
}
