#' Read tracking and feature matrices from the Caltech Fly Tracker
#'
#' @param f Path to a matlab .mat file
#' @seealso
#' \url{http://www.vision.caltech.edu/Tools/FlyTracker/documentation.html}
#' @export
#' @return For \code{read_caltech_feat} an array with dimensions time x features
#'   x flies
#' @examples
#' \dontrun{
#' feats=read_caltech_feat('~/projects/Dana/tracking/samples/210517_006-feat.mat')
#' plot(ts(feats[,"vel", "fly001"], start=0, frequency=30))
#' }
read_caltech_feat <- function(f) {
  x=rmatio::read.mat(path.expand(f))
  checkmate::assert_list(x, len=1)
  feat.names=unlist(x$feat[[1]], use.names = F)
  checkmate::assert_character(feat.names)

  feat.units=unlist(x$feat$units, recursive = F)
  feat.units <- sapply(feat.units, function(y) {
    uy=unlist(y, use.names = F)
    ifelse(length(uy)==0, NA_character_ ,as.character(uy))
  })

  checkmate::assert_character(feat.units, len = length(feat.names))
  feat.data=x$feat$data[[1]]
  # [n_flies x n_frames x n_fields double] feature data matrix
  checkmate::assert_array(feat.data, d=3L)
  nflies=dim(feat.data)[1]
  ntimepoints=dim(feat.data)[2]
  nfeats=dim(feat.data)[3]
  stopifnot(all.equal(nfeats, length(feat.names)))
  dimnames(feat.data)=list(sprintf("fly%03d", seq_len(nflies)),
                           NULL, feat.names)
  # permute from fly, time, features => time, features, fly
  feat.data=aperm(feat.data, c(2,3,1))
  attr(feat.data, "feat.names") <- feat.names
  attr(feat.data, "feat.units") <- feat.units
  feat.data
}


#' @return For \code{read_caltech_track} an array with dimensions n_frames x
#'   n_fields x n_flies and attributes:
#'
#'   \itemize{
#'
#'   \item \code{track.mat} raw tracking data (e.g. position, orientation, left
#'   wing angle)
#'
#'   \item \code{trk.names} 1 x n_fields cell names of fields in trk.data
#'
#'   \item \code{trk.flags} n_flags x 6 potential identity swaps (fly1 fly2
#'   start_fr #' end_fr ambig) The 6th column is mysterious.
#'
#'   }
#' @export
#' @rdname read_caltech_feat
#' @examples
#' \dontrun{
#' track=read_caltech_track('~/projects/Dana/tracking/samples/210517_006-track.mat')
#' plot(ts(track[,c("pos x", "pos y"), "fly001"], start=0, frequency=30))
#'
#' # plot 2 flies
#' flies=attr(track, 'flies_in_chamber')
#' plot(track[,c("pos x", "pos y"), flies[[1]][1]], type='l', col='red', asp=1)
#' lines(track[,c("pos x", "pos y"), flies[[1]][2]], type='l', col='blue')
#' }
read_caltech_track <- function(f) {
  x=rmatio::read.mat(path.expand(f))
  checkmate::assert_list(x, len=1)
  field.names=unlist(x$trk[[1]], use.names = F)
  checkmate::assert_character(field.names, min.len = 1)

  trk.flags=x$trk$flags[[1]]
  checkmate::assert_matrix(trk.flags, ncols = 6)
  colnames(trk.flags)=c("fly1", "fly2", "start_fr", "end_fr", "ambig", "X")
  trk.flags=as.data.frame(trk.flags)
  for(i in c(1:4, 6L))
    mode(trk.flags[[i]])='integer'

  trk.data=x$trk$data[[1]]
  checkmate::assert_array(trk.data, d=3L)
  nflies=dim(trk.data)[1]
  dimnames(trk.data)=list(sprintf("fly%03d", seq_len(nflies)),
                           NULL, field.names)
  # n_flies x n_frames x n_fields => n_frames x n_fields x n_flies
  trk.data=aperm(trk.data, c(2,3,1))
  attr(trk.data, 'trk.names')=field.names
  attr(trk.data, 'trk.flags')=trk.flags
  if(!is.null(x$trk$flies_in_chamber))
    attr(trk.data, 'flies_in_chamber')=unlist(x$trk$flies_in_chamber, recursive = F)
  trk.data
}
