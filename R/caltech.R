#' Read the feature matrix from the Caltech Fly Tracker
#'
#' @param f Path to a matlab .mat file
#' @seealso
#' \url{http://www.vision.caltech.edu/Tools/FlyTracker/documentation.html}
#' @export
#' @examples
#' \dontrun{
#' x=read_caltech_feat('~/projects/Dana/tracking/samples/210517_006-feat.mat')
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
