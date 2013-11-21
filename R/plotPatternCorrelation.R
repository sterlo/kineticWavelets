setMethod("plotPatternCorrelation","list",
    function(baseCorrelation){
    if(is.null(DNAPattern)){
        stop("DNAPattern is null cannot plotPatternCorrelation")
    } 
	interp=unlist(baseCorrelation$interp)[1:2^maxpwr]
    baseCorrelation=cbind(baseCorrelation$baseCorrelation[,1],baseCorrelation$baseCorrelation[,2],interp)
    names=c("IPD",'Inserts',
    paste('Pattern (',round(sum(interp)/length(interp),2),")",sep=""))
    colnames(baseCorrelation)=names
    par(mfcol=c(length(names),length(names)),mar=c(0,0,0,0.5),oma=c(4,3,3,2),xaxt="s",cex=0.5,las=1)

    levels=11

    vals=2^(1:levels)

    for(i in 1:length(names)){
		for(j in 1:length(names)){
		    temp <- vector(length=length(levels));

        if(i==j) {
            wave1=wd(baseCorrelation[,i],filter.number=1,family="DaubExPhase")
            for ( level in 1:levels){
                temp[level]=sum((accessD(wave1,level=(wave1$nlevels - level)))^2);
            }
            temp <- temp/sum(temp);
            plot(temp,type="o",col="black",ylim=c(0,0.5),axes=F);
                box(col="black",lwd=1.5);
        }
        
        if( i < j){
            
            wave1=wd(baseCorrelation[,i],filter.number=1,family="DaubExPhase")
            wave2=wd(baseCorrelation[,j],filter.number=1,family="DaubExPhase")
            for ( level in 1:levels){
                test <- cor.test((accessD(wave1,level=(wave1$nlevels - level))),                          		(accessD(wave2,level=(wave2$nlevels - level))));
                temp[level]=test$estimate
            }

            plot(temp,type="l",ylim=c(-1,1),axes=F,xlab="",ylab="");
               abline(0,0,col="grey",lty=2);
            box(col="grey"); grid(ny=0,nx=NULL);

            for(level in 1:levels){
                if(temp[level]>0) points(level,temp[level],ty='o',col="red")
                if(temp[level]<0) points(level,temp[level],ty='o',col="blue")
            }
        }
        #smooth top right
        if( i > j){
            wave1=wd(baseCorrelation[,i],filter.number=1,family="DaubExPhase")
            wave2=wd(baseCorrelation[,j],filter.number=1,family="DaubExPhase")
            for ( level in 1:levels){
                test <- cor.test((accessC(wave1,level=(wave1$nlevels - level))),                          		(accessC(wave2,level=(wave2$nlevels - level))));
                temp[level]=test$estimate
            }
            plot(temp,type="l",ylim=c(-1,1),axes=F,xlab="",ylab="");
                   abline(0,0,col="grey",lty=2);
                   box(col="grey"); grid(ny=0,nx=NULL);
            for(level in 1:levels){
                if(temp[level]>0) points(level,temp[level],ty='o',col="red")
                if(temp[level]<0) points(level,temp[level],ty='o',col="blue")
            }
        }
        labels = as.character(2^(1:levels));
        if (i==1) mtext(names[j],side=2,line=0.5,las=3,cex=1.3);
        if (j==1) mtext(names[i],side=3,line=0.5,cex=1.3);
        if (j==length(names)) axis(1,at=1:levels,labels=labels,las=3);
        if (i==length(names) & j <length(names))axis(4,tck=0.05,cex=0.5) 
        }
    }
})

