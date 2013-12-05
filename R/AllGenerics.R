#' baseCorrelation generic
#'
#' Reads a KineticWavelets object and perform base correlation
#' @param KineticWavelets \code{linkS4class{KineticWavelets}} object
#' @rdname baseCorrelation-methods
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
#' @return A list containining the kinetics for matched sequences

setGeneric("wave128Window",
           function(KineticWavelets,DNAPattern,...){
            standardGeneric("wave128Window")})

#' plotBaseCorrelation
#'
#' Reads a BaseCorrelation object and creates a plot showing
#' correlation between the smooth (top right) and detail (bottom left)
#' wavelet coeffcients at each scale, between each pair of variables.
#' The relationships between inter-pulse duration, the number of inserts,
#' and each individual nucleotide base are all examined.
#'
#' @param BaseCorrelation \code{linkS4class{BaseCorrelation}} object
#' @rdname plotBaseCorrelation-methods
#'
setGeneric('plotBaseCorrelation',
            function(BaseCorrelation,...){
            standardGeneric('plotBaseCorrelation')})
#' plotPatternCorrelation
#'
#' Reads a BaseCorrelation object and creates a plot showing
#' correlation of kinetics with aa DNA pattern. Returns an error 
#' if baseCorrelation was called without a DNAPattern
#'
#' @param BaseCorrelation \code{linkS4class{BaseCorrelation}} object
#' @rdname plotPatternCorrelation-methods
#'
setGeneric('plotPatternCorrelation',
            function(BaseCorrelation,...){
            standardGeneric('plotPatternCorrelation')})
#' plotSmoothAverage 
#' 
#' Reads a Wave128 object and creates a plot using
#' the smooth values.
#' 
#' @param Wave128 \code{linkS4class{Wave128}} object
#' @rdname plotSmoothAverage-methods
setGeneric('plotSmoothAverage',
            function(Wave128,...){
            standardGeneric('plotSmoothAverage')})
#' plotDetailAverage
#' 
#' Reads a Wave128 object a creates a plot using the 
#' detail values
#' 
#' @param Wave128 \code{linkS4class{Wave128}} object
#' @rdname plotSmoothAverage-methods

setGeneric('plotDetailAverage',
            function(Wave128,...){
            standardGeneric('plotDetailAverage')})
