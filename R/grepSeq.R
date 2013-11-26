grepSeq <- function(seqz,patternz,return.counts=T){
    seq.grep <- gregexpr(patternz,seqz)
    grep.pos=list()
    if( seq.grep[[1]][1]!=-1){
        for(i in 1:length(seq.grep[[1]])) {
            grep.pos[[i]] <- c(as.numeric(seq.grep[[1]][i]):c(-1+as.numeric(seq.grep[[1]][i])+as.numeric(attr(seq.grep[[1]], "match.length")[[i]])))
    }

    grep.seq <- c(rep(0,nchar(seqz)))
    grep.seq[unlist(grep.pos)] <- 1
    if(return.counts==T){cat("\n",length(seq.grep[[1]]),"regions found \n")}
    return(grep.seq)
    }
    if( seq.grep[[1]][1]==-1 ) {cat("\npattern not found \n")} 
}
