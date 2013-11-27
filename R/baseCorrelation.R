setMethod("baseCorrelation","KineticWavelets"
        ,
		function(KineticWavelets,DNAPattern=NULL,minReadLength=200,maxReads=10000){
    interp = NULL
    h5 = KineticWavelets@h5
    reff =KineticWavelets@reff
	if(!is.null(DNAPattern)){

        tempst <- getTemplateStrand(h5)

        endz <- getTemplateEnd(h5)
        startz <- getTemplateStart(h5)

        wFwd <- which(tempst=="0")
        wRev <- which(tempst=="1")

        fwdStart <- startz[wFwd]
        revStart <- startz[wRev]

        fwdEnd <- endz[wFwd]
        revEnd <- endz[wRev]
    ##Make reference information, assumes reff in fasta	
        reffFwd <- reff[[1]]
        reffFwd <- DNAString(c2s(reffFwd))
        reffRev <- reverseComplement(reffFwd)

        #grep for positions
        grepF <- gregexpr(DNAPattern,reffFwd)
        grepR <- gregexpr(DNAPattern,reffRev)

        gF128=list()
        gR128=list()

        if (grepF[[1]][1]!=-1){
            for(i in 1:length(grepF[[1]])) {	
                gF128[[i]] <-  c(grepF[[1]][i]+1,
                grepF[[1]][i] + attr(grepF[[1]],"match.length")[[i]])-1
            }
        }

        if (grepR[[1]][1]!=-1){
            for(i in 1:length(grepR[[1]])) {
                gR128[[i]] <- c( nchar(reffRev) - grepR[[1]][i] + 1 ,
                nchar(reffRev) - grepR[[1]][i]-attr(grepR[[1]],"match.length")[[i]] +2 
             )
        }
        }

        cat("\n checking pattern in forward sequence \n")
        interp.1 <- grepSeq(reffFwd,DNAPattern)
        cat("\n checking pattern in reverse sequence \n")
        interp.0 <- rev(grepSeq(reffRev,DNAPattern))
        totalElements=c(length(grepF[[1]][grepF[[1]]!=-1])+length(grepR[[1]][grepR[[1]]!=-1]))

        meanElementSize=sum(
        c(
        if(attr(grepF[[1]],"match.length")[1]!=-1) attr(grepF[[1]],"match.length")
        ,
        if(attr(grepR[[1]],"match.length")[1]!=-1) attr(grepR[[1]],"match.length"))
        )/totalElements


        cat("\n Reference contains",totalElements,"unique matches",sep=" ")

        rlz<-getReadLength(h5)

        fwd.rlz=rlz[wFwd]
        rev.rlz=rlz[wRev]

        g0read=list()

        if (grepR[[1]][1]!=-1){
            for(j in 1:length(gR128)){

                g0read[[j]] <- wFwd[which(fwdStart<=gR128[[j]][2] & fwdEnd>=gR128[[j]][1] & fwd.rlz>minReadLength)]
            }
        }
        ##Find reads that have windowReverse

        g1read=list()

        if (grepF[[1]][1]!=-1){
            for(j in 1:length(gF128)){

                    g1read[[j]] <- wRev[which(revStart<=gF128[[j]][1] & revEnd>=gF128[[j]][2] & rev.rlz>minReadLength) ]
            }
        }

        readsU=unique(c(unlist(g0read),unlist(g1read)))
        readsUsed=length(readsU)
        cat("\n \n File contains",readsUsed,"reads covering pattern \n",sep=" ")

        if (maxReads<readsUsed) {cat("\n Only",maxReads,"reads used.")}

            idxx=sample(unique(c(unlist(g0read),unlist(g1read))),min(maxReads,readsUsed))
        }
    if(is.null(DNAPattern)){
        rlz<-getReadLength(h5)
        whichLong=which(rlz>minReadLength)
        idxx=sample(whichLong,min(maxReads,length(whichLong)))
    }
    ipd <- getIPD (h5, idx=idxx)
    align <- getAlignments(h5 , idx=idxx)
    instsCount=list()
    if(!is.null(DNAPattern)) interp=list();
        for ( i in 1:length(ipd)){
            instsCount[[i]]=c(rep(0,length(ipd[[i]])))
            ins <- which(align[[i]][,2]=="-")
            if(length(ins)==0) {
                align[[i]]=align[[i]][,2]
            }
            if(length(ins)>0){	
                tempICount=1
                for ( j in 1:length(ins)){
                    pos=ins[j]
                    nxtup=pos+1
                        if(ins[j+1]>nxtup & j<length(ins)){
                            instsCount[[i]][nxtup]=tempICount
                            tempICount=1
                        }
                        else{
                            tempICount=tempICount+1
                        }
                }
                
                align[[i]]=align[[i]][-ins,2]
                ipd[[i]]=ipd[[i]][-ins]	
                instsCount[[i]]=instsCount[[i]][-ins]
                            
            }

            if(!is.null(DNAPattern)){
            interp[[i]]=rev(grepSeq(reverseComplement(DNAString(c2s(align[[i]]))), DNAPattern,return.counts=F))
            }

        }

    ipd=unlist(ipd)
    ipd[is.na(ipd)]=0
    align=unlist(align)
    IC=unlist(instsCount)
    maxpwr=floor(log(length(unlist(ipd)),base=2))
    ipd=ipd[1:(2^maxpwr)]
    align=align[1:(2^maxpwr)]
    IC=IC[1:(2^maxpwr)]
    aa=rep(0,2^maxpwr)
    tt=rep(0,2^maxpwr)
    cc=rep(0,2^maxpwr)
    gg=rep(0,2^maxpwr)
    aa[align=="T"]=1
    tt[align=="A"]=1
    cc[align=="G"]=1
    gg[align=="C"]=1

    tot.length=length(ipd)

    p.a=sum(aa)/tot.length
    p.t=sum(tt)/tot.length
    p.c=sum(cc)/tot.length
    p.g=sum(gg)/tot.length

    baseCorr=cbind(ipd,IC,aa,tt,cc,gg)
    names=c("IPD",'Inserts',
    paste("A (",round(p.a,2),")",sep=""),
    paste("T (",round(p.t,2),")",sep=""),
    paste("C (",round(p.c,2),")",sep=""),
    paste('G (',round(p.g,2),")",sep=""))

    colnames(baseCorr)=names

    # Maybe move this function
    return(list(baseCorrelation=baseCorr,DNAPattern=DNAPattern,interp=interp))
})
#pdf('base_corr.pdf')
#
#par(mfcol=c(length(names),length(names)),mar=c(0,0,0,0.5),oma=c(4,3,3,2),xaxt="s",cex=0.5,las=1)
#
#levels=11
#
#vals=2^(1:levels)
#
#for(i in 1:length(names)){
#		for(j in 1:length(names)){
#
#		temp <- vector(length=length(levels));
#
#
##power
#if(i==j) {
#	
#	wave1=wd(baseCorr[,i],filter.number=1,family="DaubExPhase")
#	
#	for ( level in 1:levels){
#	temp[level]=sum((accessD(wave1,level=(wave1$nlevels - level)))^2);
#	}
#        		temp <- temp/sum(temp);
#     			plot(temp,type="o",col="black",ylim=c(0,0.5),axes=F);
#         			box(col="black",lwd=1.5);
#	#axis(side=2)
#	}
#	
#	#detail bottom left
#	if( i < j){
#		
#	wave1=wd(baseCorr[,i],filter.number=1,family="DaubExPhase")
#	wave2=wd(baseCorr[,j],filter.number=1,family="DaubExPhase")
#	for ( level in 1:levels){
#							test <- cor.test((accessD(wave1,level=(wave1$nlevels - level))),                          		(accessD(wave2,level=(wave2$nlevels - level))));
#
#temp[level]=test$estimate
#}
#
#plot(temp,type="l",ylim=c(-1,1),axes=F,xlab="",ylab="");
#       abline(0,0,col="grey",lty=2);
#       box(col="grey"); grid(ny=0,nx=NULL);
#
#for(level in 1:levels){
#	if(temp[level]>0) points(level,temp[level],ty='o',col="red")
#	if(temp[level]<0) points(level,temp[level],ty='o',col="blue")
#}
#
#}
#
#	#smooth top right
#	if( i > j){
#		
#			
#	wave1=wd(baseCorr[,i],filter.number=1,family="DaubExPhase")
#	wave2=wd(baseCorr[,j],filter.number=1,family="DaubExPhase")
#	for ( level in 1:levels){
#							test <- cor.test((accessC(wave1,level=(wave1$nlevels - level))),                          		(accessC(wave2,level=(wave2$nlevels - level))));
#
#temp[level]=test$estimate
#}
#
#plot(temp,type="l",ylim=c(-1,1),axes=F,xlab="",ylab="");
#       abline(0,0,col="grey",lty=2);
#       box(col="grey"); grid(ny=0,nx=NULL);
#
#for(level in 1:levels){
#	if(temp[level]>0) points(level,temp[level],ty='o',col="red")
#	if(temp[level]<0) points(level,temp[level],ty='o',col="blue")
#}
#		
#	}
#
#
#
#
#labels = as.character(2^(1:levels));
#    			if (i==1) mtext(names[j],side=2,line=0.5,las=3,cex=1.3);
#    			if (j==1) mtext(names[i],side=3,line=0.5,cex=1.3);
#    			if (j==length(names)) axis(1,at=1:levels,labels=labels,las=3);
#	if (i==length(names) & j==1) axis(4,tck=0.05,cex=0.9) 
#	if( i != j ) axis(4,tck=0.05,cex=0.5,labels=FALSE)
#	}
#}
#
#print(length(rlz))
#print(maxpwr)
#
#dev.off()
#
#
#if(!is.null(DNAPattern)){ 
#	
#	interp=unlist(interp)[1:2^maxpwr]
#
#
#
#baseCorr=cbind(ipd,IC,interp)
#names=c("IPD",'Inserts',
#paste('Pattern (',round(sum(interp)/length(interp),2),")",sep=""))
#
#colnames(baseCorr)=names
#
#
#pdf('pattern_corr.pdf')
#
#par(mfcol=c(length(names),length(names)),mar=c(0,0,0,0.5),oma=c(4,3,3,2),xaxt="s",cex=0.5,las=1)
#
#levels=11
#
#vals=2^(1:levels)
#
#for(i in 1:length(names)){
#		for(j in 1:length(names)){
#
#		temp <- vector(length=length(levels));
#
#
##power
#if(i==j) {
#	
#	wave1=wd(baseCorr[,i],filter.number=1,family="DaubExPhase")
#	
#	for ( level in 1:levels){
#	temp[level]=sum((accessD(wave1,level=(wave1$nlevels - level)))^2);
#	}
#        		temp <- temp/sum(temp);
#     			plot(temp,type="o",col="black",ylim=c(0,0.5),axes=F);
#         			box(col="black",lwd=1.5);
#	#axis(side=2)
#	}
#	
#	#detail bottom left
#	if( i < j){
#		
#	wave1=wd(baseCorr[,i],filter.number=1,family="DaubExPhase")
#	wave2=wd(baseCorr[,j],filter.number=1,family="DaubExPhase")
#	for ( level in 1:levels){
#							test <- cor.test((accessD(wave1,level=(wave1$nlevels - level))),                          		(accessD(wave2,level=(wave2$nlevels - level))));
#
#temp[level]=test$estimate
#}
#
#plot(temp,type="l",ylim=c(-1,1),axes=F,xlab="",ylab="");
#       abline(0,0,col="grey",lty=2);
#       box(col="grey"); grid(ny=0,nx=NULL);
#
#for(level in 1:levels){
#	if(temp[level]>0) points(level,temp[level],ty='o',col="red")
#	if(temp[level]<0) points(level,temp[level],ty='o',col="blue")
#}
#
#}
#
#	#smooth top right
#	if( i > j){
#		
#			
#	wave1=wd(baseCorr[,i],filter.number=1,family="DaubExPhase")
#	wave2=wd(baseCorr[,j],filter.number=1,family="DaubExPhase")
#	for ( level in 1:levels){
#							test <- cor.test((accessC(wave1,level=(wave1$nlevels - level))),                          		(accessC(wave2,level=(wave2$nlevels - level))));
#
#temp[level]=test$estimate
#}
#
#plot(temp,type="l",ylim=c(-1,1),axes=F,xlab="",ylab="");
#       abline(0,0,col="grey",lty=2);
#       box(col="grey"); grid(ny=0,nx=NULL);
#
#for(level in 1:levels){
#	if(temp[level]>0) points(level,temp[level],ty='o',col="red")
#	if(temp[level]<0) points(level,temp[level],ty='o',col="blue")
#}
#		
#	}
#
#
#
#
#labels = as.character(2^(1:levels));
#    			if (i==1) mtext(names[j],side=2,line=0.5,las=3,cex=1.3);
#    			if (j==1) mtext(names[i],side=3,line=0.5,cex=1.3);
#    			if (j==length(names)) axis(1,at=1:levels,labels=labels,las=3);
#	if (i==length(names) & j <length(names))axis(4,tck=0.05,cex=0.5) 
#	}
#}
#
#
#dev.off()
#
#
#}
#
#})
