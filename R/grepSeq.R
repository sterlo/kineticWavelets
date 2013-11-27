grepSeq <- function(seqz,patternz,returnCounts=T){
    seqGrep <- gregexpr(patternz,seqz)
    wave.pos=list()
    if( seqGrep[[1]][1]!=-1){
        for(i in 1:length(seqGrep[[1]])) {
            wave.pos[[i]] <- c(as.numeric(seqGrep[[1]][i]):c(-1+as.numeric(seqGrep[[1]][i])+as.numeric(attr(seqGrep[[1]], "match.length")[[i]])))
    }

    wave.seq <- c(rep(0,nchar(seqz)))
    wave.seq[unlist(wave.pos)] <- 1
    if(returnCounts==T){cat("\n",length(seqGrep[[1]]),"regions found \n")}
    return(wave.seq)
    }
    if( seqGrep[[1]][1]==-1 ) {cat("\npattern not found \n")} 
}
