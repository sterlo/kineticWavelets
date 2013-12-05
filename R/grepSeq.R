#' Searches sequences for a given pattern.
#' 
#' The function returns a vector the length of the original sequence (seqz)
#' that contains a 1 if the position matched the pattern
#' otherwise 0.
#'
#' @param seqz Character vector representing the to search for a pattern.
#' @param patternz Regular expression to search for
#' @param returnCounts if \code{TRUE} print the number of matches.
#' 
#' @return A vector containing the matches.
#' @author Sterling Sawaya \email{sterlingsawaya@@gmail.com}
#' @export
#' @examples
#' grepSeq('ATCTCTCTTCTCTTC','TC')
#' # 6 regions found 
#' # [1] 0 1 1 1 1 1 1 0 1 1 1 1 0 1 1
#'

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
