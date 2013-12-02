#' Constructor for KineticWavelets class
#' @param h5 PacBio h5 Kinetics data file
#' @param reff Reference fasta file
#' @return new KineticWavelets \code{linkS4class{KineticWavelets}}
#' @author Sterling Sawaya \email{sterlingsawaya@@gmail.com}
#' @export
KineticWavelets  <- function(h5,reff){
	h5=PacBioCmpH5(h5)
	reff=read.fasta(reff)
	return(new("KineticWavelets",h5=h5,reff=reff))
}
