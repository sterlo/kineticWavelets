<<<<<<< HEAD
grepSeq <- function(seqz,patternz,returnCounts=T){
    seqGrep <- gregexpr(patternz,seqz)
    wavePos=list()
    if( seqGrep[[1]][1]!=-1){
        for(i in 1:length(seqGrep[[1]])) {
            wavePos[[i]] <- c(as.numeric(seqGrep[[1]][i]):c(-1+as.numeric(seqGrep[[1]][i])+as.numeric(attr(seqGrep[[1]], "match.length")[[i]])))
    }

    waveSeq <- c(rep(0,nchar(seqz)))
    waveSeq[unlist(wavePos)] <- 1
    if(returnCounts==T){cat("\n",length(seqGrep[[1]]),"regions found \n")}
    return(waveSeq)
    }
    if( seqGrep[[1]][1]==-1 ) {cat("\npattern not found \n")} 
}
