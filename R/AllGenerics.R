#' baseCorrelation generic
#'
#' Reads a KineticWavelets object and perform base correlation
#' @param KineticWavelets \code{linkS4class{KineticWavelets}} object
#' @rdname baseCorrelation-methods
#' @export
#' @author Sterling Sawaya \email{sterlingsawaya@@gmail.com}
#' @return A list containining the base correlations

setGeneric("baseCorrelation",
           function(KineticWavelets,...){standardGeneric("baseCorrelation")})

#' wave128Window generic
#'
#' Reads a KineticWavelets object and creates wavelets     for reads that match the DNAPattern
#' @param KineticWavelets \code{linkS4class{KineticWavelets}} object
#' @param DNAPattern The DNA pattern to search for in the reads
#' @rdname wave128Window-methods
#' @author Sterling Sawaya \email{sterlingsawaya@@gmail.com}
#' @export
#' @return A list containining the kinetics for matched sequences

setGeneric("wave128Window",
           function(KineticWavelets,DNAPattern,...){standardGeneric("wave128Window")})
