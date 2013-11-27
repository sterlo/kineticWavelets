KineticWavelets  <- function(h5,reff){
	h5=PacBioCmpH5(h5)
	reff=read.fasta(reff)
	return(new("KineticWavelets",h5=h5,reff=reff))
};