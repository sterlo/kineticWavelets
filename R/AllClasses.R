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
#' Wave128 Class.
#'
#' A Wave128 class is used to track the results of wave128Window function run
#' Object is to be passed to the plotting functions.
#' 
#' @section Slots:
#'     \describe{
#'              \item{\code{detailWave}:}{Object of class \code{"list"}.}
#'              \item{\code{smoothWave}:}{Object of class \code{"list"}.}
#'              \item{\code{totalElements}:}{Object of class \code{"numeric"}.}
#'              \item{\code{meanElementSize}:}{Object of class \code{"numeric"}.}
#'              \item{\code{alignments}:}{Object of class \code{"list"}.}
#'              \item{\code{DNAPattern}:}{Object of class \code{"character"}.}
#'              \item{\code{shiftWindow}:}{Object of class \code{"numeric"}.}
#'              \item{\code{shrink}:}{Object of class \code{"numeric"}.}
#'      }
#' @exportClass Wave128     
setClass('Wave128',
        representation(detailWave ='list',smoothWave='list',                      totalElements='numeric',meanElementSize='numeric',alignments='list',DNAPattern='character',shiftWindow='numeric',shrink='numeric'),
        prototype(detailWave=NULL,smoothWave=NULL,
                    totalElements=NULL,meanElementSize=NULL,alignments=NULL,DNAPattern=NULL,shiftWindow=NULL,
            shrink=NULL)
    

);
#' BaseCorrelation Class.
#'
#' A BaseCorrelation class is sude to track the results of the baseCorrelation function.
#'
#' @section Slots:
#'   \describe{
#'          \item{\code{base}:}{Object of class \code{"matrix"}.}
#'          \item{\code{DNAPattern}:}{Object of class \code{"character"}}
#'          \item{\code{interp}:}{Object of class \code{"list"}}
#'   }
#' @exportClass BaseCorrelation


setClass('BaseCorrelation',
        representation(baseCorr='matrix',DNAPattern='character',interp='list'),
        prototype(baseCorr=NULL,DNAPattern=NULL,interp=NULL)
);
