<<<<<<< HEAD
grepSeq <- function(seqz,patternz,returnCounts=T){
    seqGrep <- gregexpr(patternz,seqz)
    grepPos=list()
    if( seqGrep[[1]][1]!=-1){
        for(i in 1:length(seqGrep[[1]])) {
            grepPos[[i]] <- c(as.numeric(seqGrep[[1]][i]):c(-1+as.numeric(seqGrep[[1]][i])+as.numeric(attr(seqGrep[[1]], "match.length")[[i]])))
    }

    grepSeq <- c(rep(0,nchar(seqz)))
    grepSeq[unlist(grepPos)] <- 1
    if(returnCounts==T){cat("\n",length(seqGrep[[1]]),"regions found \n")}
    return(grepSeq)
    }
    if( seqGrep[[1]][1]==-1 ) {cat("\npattern not found \n")} 
}
