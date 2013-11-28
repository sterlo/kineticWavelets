<<<<<<< HEAD
#' Reads a KineticWavelets object and creates wavelets for reads that match the DNAPattern
#' @param KineticWavelets \code{linkS4class{KineticWavelets}} object
#' @param DNAPattern The DNA pattern to search for in the reads.
#' @param maxReads The maximum number of reads for each region
#' @param shiftWindow 
#' @param totalTime 
#' @param filterNumber 
#' @param shrink
#' @return A list containing all the wavelet correlations to be passed to plotting 
#'         functions or investigated.
#' @author Sterling Sawaya \email{sterlingsawaya@gmail.com}
#' @examples
#' kin = KineticWavelets(h5,reff)
#' wav = waveCorrelation(kin,DNAPattern="GCGCGCGCG")
#'

setMethod("waveCorrelation","KineticWavelets",
            function(KineticWavelets,DNAPattern,maxReads=1000,shiftWindow=64,totalTime=TRUE,filterNumber=1,shrink=1){
    h5 = KineticWavelets@h5
    reff = KineticWavelets@reff

    tempst <- getTemplateStrand(h5)

    endz <- getTemplateEnd(h5)
    startz <- getTemplateStart(h5)

    w.fwd <- which(tempst=="0")
    w.rev <- which(tempst=="1")

    fwd.start <- startz[w.fwd]
    rev.start <- startz[w.rev]

    fwd.end <- endz[w.fwd]
    rev.end <- endz[w.rev]

    reff.fwd <- reff[[1]]
    reff.fwd <- DNAString(c2s(reff.fwd))
    reff.rev <- reverseComplement(reff.fwd)

    grepF <- gregexpr(DNAPattern,reff.fwd)
    grepR <- gregexpr(DNAPattern,reff.rev)

    gF128=list()
    gR128=list()
    forward_and_reverse_reads = matrix(nrow=length(grepF[[1]])+ length(grepR[[1]]),ncol=2)
    no_forward = length(grepF[[1]])
    no_rev = length(grepR[[1]])

    read_number = 1
    if (grepF[[1]][1]!=-1){
    for(i in 1:length(grepF[[1]])) {
        temp <- c(grepF[[1]][i] + attr(grepF[[1]],"match.length")[[i]])-1
        gF128[[i]] <- c(c(temp - 127 + shiftWindow),c(temp + shiftWindow))
        forward_and_reverse_reads[read_number,] = c(c(temp - 127 + shiftWindow),c(temp + shiftWindow))
        read_number = read_number + 1
        }
    }

    if (grepR[[1]][1]!=-1){
    for(i in 1:length(grepR[[1]])) {
        
        temp <- nchar(reff.rev) - grepR[[1]][i]-attr(grepR[[1]],"match.length")[[i]] + 2
        gR128[[i]] <- c(temp-shiftWindow,temp+127-shiftWindow)
        forward_and_reverse_reads[read_number,] = c(temp-shiftWindow,temp+127-shiftWindow)
        read_number = read_number + 1
    }
    }
    cat("\n checking pattern in forward sequence \n")
    interp.1 <- grepSeq(reff.fwd,DNAPattern)
    cat("\n checking pattern in reverse sequence \n")
    interp.0 <- rev(grepSeq(reff.rev,DNAPattern))


    totalElements=c(length(grepF[[1]][grepF[[1]]!=-1])+length(grepR[[1]][grepR[[1]]!=-1]))

    meanElementSize=sum(
    c(
    if(attr(grepF[[1]],"match.length")[1]!=-1) attr(grepF[[1]],"match.length")
    ,
    if(attr(grepR[[1]],"match.length")[1]!=-1) attr(grepR[[1]],"match.length"))
    )/totalElements


    cat("\n Reference contains",totalElements,"unique matches with potentially overlapping 128 bp windows.  ",sep=" ")

    ##Find reads that have window forward



    gread=list()

    if (grepR[[1]][1]!=-1){
    for(j in 1:length(gR128)){

            gread[[j]] <- w.fwd[which(fwd.start<=gR128[[j]][1] & fwd.end>=gR128[[j]][2])]
        }
    }
    ##Find reads that have window reverse
    fwd_length = length(gread)

    if (grepF[[1]][1]!=-1){
    for(k in 1:length(gF128)){

            gread[[k+length(gR128)]] <- w.rev[which(rev.start<=gF128[[k]][1] & rev.end>=gF128[[k]][2])]
        }
    }
    readsUsed=c(length(unlist(gread)))
    cat("\n \n File contains",readsUsed,"reads that cover 128 bp windows of interest.\n",sep=" ")

    if (maxReads<readsUsed) {cat("\n Only the first",maxReads,"reads for each region are used. \n")}

    smoothWave=list()
    k = 1
    ipd_nums = vector(length=length(gread),mode="integer")
    for ( i in 1:length(gread)){
        idx = gread[[i]]
        if(length(idx) > 0){
        ipd = getIPD(h5,idx=idx)
        ipd_nums[i] = min(length(ipd),maxReads) 
        }
        else{
            ipd_nums[i] = 0
        }
    }
    no_reads = sum(ipd_nums) * 128
    for ( j in 1:8){
        smoothWave[[j]] = matrix(nrow=(no_reads), ncol=2)
    }
    detailWave=list()
    for ( j in 1:3){
        detailWave[[j]] = matrix(nrow=(no_reads), ncol=2)
        }

    k=1
    reading_g0 = TRUE
    read_insert = 1
    while ( k <= length(gread)){
        idx = gread[[k]]
        if(k <= fwd_length){
            reading_g0 = T
        }else{
            reading_g0 = F  
        }
        if(length(idx)>0){
        posit <- getTemplatePosition(h5,idx=idx)
        align <- getAlignments(h5,idx=idx)
        ipd <- getIPD(h5,idx=idx)
        pw <- getPulseWidth(h5,idx=idx)
        
        for ( i in 1:min(length(ipd),maxReads)){
            ins <- which(align[[i]][,2]=="-")
            
      if(length(ins)==0) {
        align[[i]]=align[[i]][,2]

    if( totalTime==TRUE) {
        instsTime=c(rep(0,length(ipd[[i]])));
        }
    }else{
        if(totalTime==TRUE){
        instsTime=c(rep(0,length(ipd[[i]])))
        tempITime=0
        for ( j in 1:length(ins)){
            pos=ins[j]
              nxtup=pos+1
            
              if(ins[j+1]>nxtup & j<length(ins)){
                  instsTime[nxtup] <- tempITime+ipd[[i]][pos]+pw[[i]][pos]
                
                  tempITime=0
              }
              else{
                  tempITime=tempITime+ipd[[i]][pos]+pw[[i]][pos]
              }
          }
        }
    ##remove insert from all data       
            
            align[[i]] <- align[[i]][-ins,2]
            ipd[[i]] <- ipd[[i]][-ins]
            posit[[i]] <- posit[[i]][-ins]
            pw[[i]] <- pw[[i]][-ins]
            if (totalTime==TRUE) instsTime <- instsTime[-ins];
            }
        if(reading_g0){
          positz <- c(gR128[[k]][1]:gR128[[k]][2])  
        } else {
            positz <- c(gF128[[k-fwd_length]][2]:gF128[[k-fwd_length]][1])
        }
        align[[i]] <- align[[i]][match(positz,posit[[i]])]
        ipd[[i]] <- ipd[[i]][match(positz,posit[[i]])]
        pw[[i]] <- pw[[i]][match(positz,posit[[i]])]
        if (totalTime==TRUE)   instsTime <- instsTime[match(positz,posit[[i]])];
        posit[[i]] <- positz
        if(reading_g0){ 
        interp <- interp.0[posit[[i]]]
        }else{
        interp <- interp.1[posit[[i]]]
        }
        ipd[[i]][is.na(ipd[[i]])] <- 0
        pw[[i]][is.na(pw[[i]])] <- 0
        
        
           if (totalTime==TRUE){
           ipd.wst <- wd(rev(ipd[[i]]+pw[[i]]+instsTime),family="DaubExPhase",filter.number=filterNumber,type="station")
       }
       else{
                ipd.wst <- wd(rev(ipd[[i]]),family="DaubExPhase",filter.number=filterNumber,type="station")
           }
           
           
    pat.wst <- wd(rev(interp),family="DaubExPhase",filter.number=filterNumber,type="station")
                    st.index = 128 * (read_insert - 1) + 1
                    end.index = 128 * (read_insert)

                                    for ( j in 1:8){
                                smoothWave[[j]][st.index:end.index,] = cbind(accessC(ipd.wst,level=j-1),accessC(pat.wst,level=j-1))
                                    
                }
        for (j in 1:3){
                detailWave[[j]][st.index:end.index,] = matrix(accessD(ipd.wst,level=7-j),ncol=1)
                    }
                read_insert = read_insert + 1   
                }
        
        }
        k=k+1
    }   
    if (shrink<1){
            for ( j in 1:8){
            qnt  <-  quantile(smoothWave[[j]][,1],probs=1-shrink)
        smoothWave[[j]][,1]=smoothWave[[j]][,1]-qnt
        smoothWave[[j]][which(smoothWave[[j]][,1]<0),1]=0
                }
                
                for ( j in 1:3){
            qnt  <-   quantile(abs(detailWave[[j]]),probs=1-shrink)
        
        detailWave[[j]][abs(detailWave[[j]])<=qnt]=0
        detailWave[[j]][detailWave[[j]]<0] <- detailWave[[j]][detailWave[[j]]<0]+qnt
        detailWave[[j]][detailWave[[j]]>0] <- detailWave[[j]][detailWave[[j]]>0]-qnt
                }   
    }
        return(list(detailWave=detailWave,smoothWave=smoothWave,totalElements=totalElements,meanElementSize=meanElementSize,alignments=align,DNAPattern=DNAPattern,shiftWindow=shiftWindow,shrink=shrink))
})
=======
#' Reads a KineticWavelets object and creates wavelets for reads that match the DNAPattern
#' @param KineticWavelets \code{linkS4class{KineticWavelets}} object
#' @param DNAPattern The DNA pattern to search for in the reads.
#' @param maxReads The maximum number of reads for each region
#' @param shiftWindow 
#' @param totalTime 
#' @param filterNumber 
#' @param shrink
#' @return A list containing all the wavelet correlations to be passed to plotting 
#'         functions or investigated.
#' @author Sterling Sawaya \email{sterlingsawaya@gmail.com}
#' @examples
#' kin = KineticWavelets(h5,reff)
#' wav = waveCorrelation(kin,DNAPattern="GCGCGCGCG")
#'

setMethod("waveCorrelation","KineticWavelets",
            function(KineticWavelets,DNAPattern,maxReads=1000,shiftWindow=64,totalTime=TRUE,filterNumber=1,shrink=1){
    h5 = KineticWavelets@h5
    reff = KineticWavelets@reff

    tempst <- getTemplateStrand(h5)

    endz <- getTemplateEnd(h5)
    startz <- getTemplateStart(h5)

    w.fwd <- which(tempst=="0")
    w.rev <- which(tempst=="1")

    fwd.start <- startz[w.fwd]
    rev.start <- startz[w.rev]

    fwd.end <- endz[w.fwd]
    rev.end <- endz[w.rev]

    reff.fwd <- reff[[1]]
    reff.fwd <- DNAString(c2s(reff.fwd))
    reff.rev <- reverseComplement(reff.fwd)

    grepF <- gregexpr(DNAPattern,reff.fwd)
    grepR <- gregexpr(DNAPattern,reff.rev)

    gF128=list()
    gR128=list()
    forward_and_reverse_reads = matrix(nrow=length(grepF[[1]])+ length(grepR[[1]]),ncol=2)
    no_forward = length(grepF[[1]])
    no_rev = length(grepR[[1]])

    read_number = 1
    if (grepF[[1]][1]!=-1){
    for(i in 1:length(grepF[[1]])) {
        temp <- c(grepF[[1]][i] + attr(grepF[[1]],"match.length")[[i]])-1
        gF128[[i]] <- c(c(temp - 127 + shiftWindow),c(temp + shiftWindow))
        forward_and_reverse_reads[read_number,] = c(c(temp - 127 + shiftWindow),c(temp + shiftWindow))
        read_number = read_number + 1
        }
    }

    if (grepR[[1]][1]!=-1){
    for(i in 1:length(grepR[[1]])) {
        
        temp <- nchar(reff.rev) - grepR[[1]][i]-attr(grepR[[1]],"match.length")[[i]] + 2
        gR128[[i]] <- c(temp-shiftWindow,temp+127-shiftWindow)
        forward_and_reverse_reads[read_number,] = c(temp-shiftWindow,temp+127-shiftWindow)
        read_number = read_number + 1
    }
    }
    cat("\n checking pattern in forward sequence \n")
    interp.1 <- grepSeq(reff.fwd,DNAPattern)
    cat("\n checking pattern in reverse sequence \n")
    interp.0 <- rev(grepSeq(reff.rev,DNAPattern))


    totalElements=c(length(grepF[[1]][grepF[[1]]!=-1])+length(grepR[[1]][grepR[[1]]!=-1]))

    meanElementSize=sum(
    c(
    if(attr(grepF[[1]],"match.length")[1]!=-1) attr(grepF[[1]],"match.length")
    ,
    if(attr(grepR[[1]],"match.length")[1]!=-1) attr(grepR[[1]],"match.length"))
    )/totalElements


    cat("\n Reference contains",totalElements,"unique matches with potentially overlapping 128 bp windows.  ",sep=" ")

    ##Find reads that have window forward



    gread=list()

    if (grepR[[1]][1]!=-1){
    for(j in 1:length(gR128)){

            gread[[j]] <- w.fwd[which(fwd.start<=gR128[[j]][1] & fwd.end>=gR128[[j]][2])]
        }
    }
    ##Find reads that have window reverse
    fwd_length = length(gread)

    if (grepF[[1]][1]!=-1){
    for(k in 1:length(gF128)){

            gread[[k+length(gR128)]] <- w.rev[which(rev.start<=gF128[[k]][1] & rev.end>=gF128[[k]][2])]
        }
    }
    readsUsed=c(length(unlist(gread)))
    cat("\n \n File contains",readsUsed,"reads that cover 128 bp windows of interest.\n",sep=" ")

    if (maxReads<readsUsed) {cat("\n Only the first",maxReads,"reads for each region are used. \n")}

    smoothWave=list()
    k = 1
    ipd_nums = vector(length=length(gread),mode="integer")
    for ( i in 1:length(gread)){
        idx = gread[[i]]
        if(length(idx) > 0){
        ipd = getIPD(h5,idx=idx)
        ipd_nums[i] = min(length(ipd),maxReads) 
        }
        else{
            ipd_nums[i] = 0
        }
    }
    no_reads = sum(ipd_nums) * 128
    for ( j in 1:8){
        smoothWave[[j]] = matrix(nrow=(no_reads), ncol=2)
    }
    detailWave=list()
    for ( j in 1:3){
        detailWave[[j]] = matrix(nrow=(no_reads), ncol=2)
        }

    k=1
    reading_g0 = TRUE
    read_insert = 1
    while ( k <= length(gread)){
        idx = gread[[k]]
        if(k <= fwd_length){
            reading_g0 = T
        }else{
            reading_g0 = F  
        }
        if(length(idx)>0){
        posit <- getTemplatePosition(h5,idx=idx)
        align <- getAlignments(h5,idx=idx)
        ipd <- getIPD(h5,idx=idx)
        pw <- getPulseWidth(h5,idx=idx)
        
        for ( i in 1:min(length(ipd),maxReads)){
            ins <- which(align[[i]][,2]=="-")
            
      if(length(ins)==0) {
        align[[i]]=align[[i]][,2]

    if( totalTime==TRUE) {
        instsTime=c(rep(0,length(ipd[[i]])));
        }
    }else{
        if(totalTime==TRUE){
        instsTime=c(rep(0,length(ipd[[i]])))
        tempITime=0
        for ( j in 1:length(ins)){
            pos=ins[j]
              nxtup=pos+1
            
              if(ins[j+1]>nxtup & j<length(ins)){
                  instsTime[nxtup] <- tempITime+ipd[[i]][pos]+pw[[i]][pos]
                
                  tempITime=0
              }
              else{
                  tempITime=tempITime+ipd[[i]][pos]+pw[[i]][pos]
              }
          }
        }
    ##remove insert from all data       
            
            align[[i]] <- align[[i]][-ins,2]
            ipd[[i]] <- ipd[[i]][-ins]
            posit[[i]] <- posit[[i]][-ins]
            pw[[i]] <- pw[[i]][-ins]
            if (totalTime==TRUE) instsTime <- instsTime[-ins];
            }
        if(reading_g0){
          positz <- c(gR128[[k]][1]:gR128[[k]][2])  
        } else {
            positz <- c(gF128[[k-fwd_length]][2]:gF128[[k-fwd_length]][1])
        }
        align[[i]] <- align[[i]][match(positz,posit[[i]])]
        ipd[[i]] <- ipd[[i]][match(positz,posit[[i]])]
        pw[[i]] <- pw[[i]][match(positz,posit[[i]])]
        if (totalTime==TRUE)   instsTime <- instsTime[match(positz,posit[[i]])];
        posit[[i]] <- positz
        if(reading_g0){ 
        interp <- interp.0[posit[[i]]]
        }else{
        interp <- interp.1[posit[[i]]]
        }
        ipd[[i]][is.na(ipd[[i]])] <- 0
        pw[[i]][is.na(pw[[i]])] <- 0
        
        
           if (totalTime==TRUE){
           ipd.wst <- wd(rev(ipd[[i]]+pw[[i]]+instsTime),family="DaubExPhase",filter.number=filterNumber,type="station")
       }
       else{
                ipd.wst <- wd(rev(ipd[[i]]),family="DaubExPhase",filter.number=filterNumber,type="station")
           }
           
           
    pat.wst <- wd(rev(interp),family="DaubExPhase",filter.number=filterNumber,type="station")
                    st.index = 128 * (read_insert - 1) + 1
                    end.index = 128 * (read_insert)

                                    for ( j in 1:8){
                                smoothWave[[j]][st.index:end.index,] = cbind(accessC(ipd.wst,level=j-1),accessC(pat.wst,level=j-1))
                                    
                }
        for (j in 1:3){
                detailWave[[j]][st.index:end.index,] = matrix(accessD(ipd.wst,level=7-j),ncol=1)
                    }
                read_insert = read_insert + 1   
                }
        
        }
        k=k+1
    }   
    if (shrink<1){
            for ( j in 1:8){
            qnt  <-  quantile(smoothWave[[j]][,1],probs=1-shrink)
        smoothWave[[j]][,1]=smoothWave[[j]][,1]-qnt
        smoothWave[[j]][which(smoothWave[[j]][,1]<0),1]=0
                }
                
                for ( j in 1:3){
            qnt  <-   quantile(abs(detailWave[[j]]),probs=1-shrink)
        
        detailWave[[j]][abs(detailWave[[j]])<=qnt]=0
        detailWave[[j]][detailWave[[j]]<0] <- detailWave[[j]][detailWave[[j]]<0]+qnt
        detailWave[[j]][detailWave[[j]]>0] <- detailWave[[j]][detailWave[[j]]>0]-qnt
                }   
    }
        return(list(detailWave=detailWave,smoothWave=smoothWave,totalElements=totalElements,meanElementSize=meanElementSize,alignments=align,DNAPattern=DNAPattern,shiftWindow=shiftWindow,shrink=shrink))
})
;
>>>>>>> bfffe25262de104313095b3e9c1b03e0448d5129
