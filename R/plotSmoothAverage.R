#' @exportMethod plotSmoothAverage
setMethod(f="plotSmoothAverage","Wave128",
function(Wave128){
        shrink = Wave128@shrink
        meanElementSize = Wave128@meanElementSize
        align = Wave128@alignments
        smoothWave = Wave128@smoothWave
        totalElements = Wave128@totalElements
        shiftWindow = Wave128@shiftWindow
        DNAPattern = Wave128@DNAPattern
par(mfcol=c(3,1))   
mean0=c()
qnt0=matrix(NA,nrow=4,ncol=128)
        for(i in 1:128){
            idx=seq(i,length(smoothWave[[8]][,1]),by=128)
            mean0[i]=mean(smoothWave[[8]][idx,1])
            qnt0[1,i]=quantile(smoothWave[[8]][idx,1],0.05)                 
            qnt0[2,i]=quantile(smoothWave[[8]][idx,1],0.95)
            qnt0[3,i]=quantile(smoothWave[[8]][idx,1],0.1)
            qnt0[4,i]=quantile(smoothWave[[8]][idx,1],0.9)
        }
        mainz=paste("IPD-wavelet for ",DNAPattern," (",nrow(smoothWave[[1]])/128,
        " reads in ",totalElements," regions)",sep="")
        poz=c()
        if(totalElements==1){poz=positz[c(1,8*(1:16))]}else{poz=c(1,8*(1:16))}
        if(shrink<1)(mainz=paste(shrink,"Shrunken mean IPD for ",DNAPattern," (",nrow(smoothWave[[8]])/128," reads/",totalElements,
        " regions)",sep=""))
        
        plot(c(1:128),rev(mean0),ty="l",bty="n",ylab="IPD (s)",xlab="position in window",main=mainz,ylim=c(0,1.1*max(qnt0,mean0)),axes=F)
        axis(side=1,labels=poz,at=c(1,8*(1:16)))
        axis(side=2)
        rect(shiftWindow+0.5,-100,shiftWindow+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
        segments(c(1:128),rev(qnt0[1,]),c(1:128),rev(qnt0[2,]),col=rgb(0,0,0,0.2))
        segments(c(1:128),rev(qnt0[3,]),c(1:128),rev(qnt0[4,]),col=rgb(0,0,0,0.4))

        if (totalElements==1){
             text(1:128,1.05*(max(qnt0,mean0)),align[[1]],cex=0.5)
        }
        if(totalElements!=1) {
            text((shiftWindow-5),1.05*(max(qnt0,mean0)),"Read:",cex=0.8)
            text((shiftWindow+1):(shiftWindow+floor(meanElementSize)),1.05*(max(qnt0,mean0)),align[[1]][(shiftWindow+1):(shiftWindow+floor(meanElementSize))],cex=0.5);
        }
            
        mean1=c()
        qnt1=matrix(NA,nrow=4,ncol=128)
        for(i in 1:128){
            idx=seq(i,length(smoothWave[[7]][,1]),by=128)
            mean1[i]=mean(smoothWave[[7]][idx,1])
            qnt1[1,i]=quantile(smoothWave[[7]][idx,1],0.05)                 
            qnt1[2,i]=quantile(smoothWave[[7]][idx,1],0.95)
            qnt1[3,i]=quantile(smoothWave[[7]][idx,1],0.1)
            qnt1[4,i]=quantile(smoothWave[[7]][idx,1],0.9)  
        }
            plot(c(1:128),rev(mean1),ty="l",bty="n",ylab=expression(E(s(2))),xlab="position in window",main="2 bp smoothing",ylim=c(0,1.1*max(qnt1,mean1)),axes=F)
            byz=c(8*(1:16))
            axis(side=1,labels=poz,at=c(1,byz))
            axis(side=2)
            
            rect(shiftWindow+0.5,-100,shiftWindow+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
            segments(c(1:128),rev(qnt1[1,]),c(1:128),rev(qnt1[2,]),col=rgb(0,0,0,0.2))
            segments(c(1:128),rev(qnt1[3,]),c(1:128),rev(qnt1[4,]),col=rgb(0,0,0,0.4))

                
        if (totalElements==1){
             text(1:128,1.05*(max(qnt1,mean1)),align[[1]],cex=0.5)
        }

        if(totalElements!=1) {
            text((shiftWindow-5),1.05*(max(qnt1,mean1)),"Read:",cex=0.8)
            text((shiftWindow+1):(shiftWindow+floor(meanElementSize)),1.05*(max(qnt1,mean1)),align[[1]][(shiftWindow+1):(shiftWindow+floor(meanElementSize))],cex=0.5);
        }
        mean2=c()
        qnt2=matrix(NA,nrow=4,ncol=128)                 
        for(i in 1:128){
            idx=seq(i,length(smoothWave[[6]][,1]),by=128)
            mean2[i]=mean(smoothWave[[6]][idx,1])
            qnt2[1,i]=quantile(smoothWave[[6]][idx,1],0.05)                 
            qnt2[2,i]=quantile(smoothWave[[6]][idx,1],0.95)
            qnt2[3,i]=quantile(smoothWave[[6]][idx,1],0.1)
            qnt2[4,i]=quantile(smoothWave[[6]][idx,1],0.9)  
        }
        
        plot(c(1:128),rev(mean2),ty="l",bty="n",ylab=expression(E(s(4))),xlab="position in window",main="4 bp smoothing",ylim=c(0,1.1*max(qnt2,mean2)),axes=F)
        axis(side=1,labels=poz,at=c(1,byz))
        axis(side=2)
        rect(shiftWindow+0.5,-100,shiftWindow+meanElementSize+.5,100,col=rgb(0,0,1,0.1),border=NA)

        segments(c(1:128),rev(qnt2[1,]),c(1:128),rev(qnt2[2,]),col=rgb(0,0,0,0.2))
        segments(c(1:128),rev(qnt2[3,]),c(1:128),rev(qnt2[4,]),col=rgb(0,0,0,0.4))
                
        if (totalElements==1){
            text(1:128,1.05*(max(qnt2,mean2)),align[[1]],cex=0.5)
        }
    
    if(totalElements!=1){
     text((shiftWindow-5),1.05*(max(qnt2,mean2)),"Read:",cex=0.8)
     text((shiftWindow+1):(shiftWindow+floor(meanElementSize)),1.05*(max(qnt2,mean2)),align[[1]][(shiftWindow+1):(shiftWindow+floor(meanElementSize))],cex=0.5)
    }
            
});
