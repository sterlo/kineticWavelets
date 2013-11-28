#' KineticWavelets class.
#' 
#' The class is used to track a CmpH5 file and a reference fasta 
#' to be used for further analysis in the KineticWavelets package.
#'
#' \section{Slots}{
#'     \describe{
#'          \item{\code{h5}:}{
#'  }
#' }
#'

setClass('KineticWavelets',
		representation(h5='PacBioCmpH5',reff='list'),
	prototype(h5=NULL,reff=NULL)
);