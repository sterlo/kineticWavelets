#' KineticWavelets class.
#' 
#' The class is used to track a CmpH5 file and a reference fasta 
#' to be used for further analysis in the KineticWavelets package.
#'
#'@section Slots:
#'     \describe{
#'          \item{\code{h5}:}{ Object of class \code{"PacBioCmpH5"}.}
#'          \item{\code{reff}:}{ Object of class \code{"list"}.}
#' }
#' @import pbh5
#' @import Biostrings
#' @import seqinr
#' @import wavethresh
#' @import ShortRead
#' @exportClass KineticWavelets

setClass('KineticWavelets',
		representation(h5='PacBioCmpH5',reff='list'),
	prototype(h5=NULL,reff=NULL)
);
