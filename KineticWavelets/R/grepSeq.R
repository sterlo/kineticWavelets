
	grepSeq <- function(seqz,patternz,return.counts=T){
		
		seq.grep <- gregexpr(patternz,seqz)

		wave.pos=list()
		if( seq.grep[[1]][1]!=-1){
		for(i in 1:length(seq.grep[[1]])) {
		wave.pos[[i]] <- c(as.numeric(seq.grep[[1]][i]):c(-1+as.numeric(seq.grep[[1]][i])+as.numeric(attr(seq.grep[[1]], "match.length")[[i]])))
	}

	wave.seq <- c(rep(0,nchar(seqz)))
	wave.seq[unlist(wave.pos)] <- 1
	if(return.counts==T){cat("\n",length(seq.grep[[1]]),"regions found \n")}
	return(wave.seq)
	}
	if( seq.grep[[1]][1]==-1 ) {cat("\npattern not found \n")} 
	}
