#######################################
#######################################
#Sterling Sawaya
#University of Otago
#March 7 2013
#######################################
#This file contains functions to make smoothed plots of wavelet coefficients for finer scales around a pattern of interest
#and to test for correlations between wavelet coefficients for any pattern

##All code has been tested to run on R 2.14.0 GUI 1.42 Leopard build 64-bit 






#h5<- PacBioCmpH5("/Users/sterling/Desktop/PacBio/Bacterial_kinetics/G.metallireducens.cmp.h5")
#reff=read.fasta("/Users/sterling/Desktop/PacBio/Bacterial_kinetics/Gmet_GS15_.fasta")
#DNAPattern="GGGGCTCCCGGCGGAGATCGGGG";shift.window=64;total.time=TRUE;max.reads=200;shrink=1;filterNumber=1;plotDetailAvg=TRUE

#DNAPattern="(GG(.){1,3}GGG(.){1,10}GGG(.){1,7}GGG)"


###################################
###################################
#Main function to calculate wavelet averages and plot
#The shift window determine where in the 128 bp window the element of interest lies.  Default is 32 so that the element starts at base 33 in the window (i.e. shifted from 1 by 32).
setGeneric("waveCorrelation",
           function(Object, ...){standardGeneric("waveCorrelation")})


setMethod("waveCorrelation","KineticWavelets",
			function(Object,DNAPattern,max.reads=1000,shifted=64){
	h5 = Object@h5
	reff = Object@reff
	DNAPattern = Object@DNAPattern
	shift.window = Object@shift.window
	total.time = Object@total.time
	max.reads = Object@max.reads
	filterNumber = Object@filterNumber
	shrink = Object@shrink
	##Function to turn sequences into numeric vectors
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

#Get cmp.h5 strands and read positions	

tempst <- getTemplateStrand(h5)

endz <- getTemplateEnd(h5)
startz <- getTemplateStart(h5)

w.fwd <- which(tempst=="0")
w.rev <- which(tempst=="1")

fwd.start <- startz[w.fwd]
rev.start <- startz[w.rev]

fwd.end <- endz[w.fwd]
rev.end <- endz[w.rev]


##Make reference information, assumes reff in fasta	
	
reff.fwd <- reff[[1]]
reff.fwd <- DNAString(c2s(reff.fwd))
reff.rev <- reverseComplement(reff.fwd)

#grep for positions
grepF <- gregexpr(DNAPattern,reff.fwd)
grepR <- gregexpr(DNAPattern,reff.rev)

##This part gets slightly complicated.  We want the read corresponding to the position of the element on the opposite strand.  This code makes the 128 bp window with a shift from the start of the element on the template strand. 

#Window position forward template
gF128=list()
#Window position reverse template
gR128=list()

# Reads forward and reverse strands
forward_and_reverse_reads = matrix(nrow=length(grepF[[1]])+ length(grepR[[1]]),ncol=2)
# length of the forward and reverser
no_forward = length(grepF[[1]])
no_rev = length(grepR[[1]])

# Set read numbers for full matrix
read_number = 1
if (grepF[[1]][1]!=-1){
for(i in 1:length(grepF[[1]])) {
	temp <- c(grepF[[1]][i] + attr(grepF[[1]],"match.length")[[i]])-1
	gF128[[i]] <- c(c(temp - 127 + shift.window),c(temp + shift.window))
	forward_and_reverse_reads[read_number,] = c(c(temp - 127 + shift.window),c(temp + shift.window))
	read_number = read_number + 1
}
}

if (grepR[[1]][1]!=-1){
for(i in 1:length(grepR[[1]])) {
	
	temp <- nchar(reff.rev) - grepR[[1]][i]-attr(grepR[[1]],"match.length")[[i]] + 2
	gR128[[i]] <- c(temp-shift.window,temp+127-shift.window)
	forward_and_reverse_reads[read_number,] = c(temp-shift.window,temp+127-shift.window)
	read_number = read_number + 1
}
}

	####################################################
#Set up interpretation of sequence

cat("\n checking pattern in forward sequence \n")
interp.1 <- waveSeq(reff.fwd,DNAPattern)
cat("\n checking pattern in reverse sequence \n")
interp.0 <- rev(waveSeq(reff.rev,DNAPattern))


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
####################################################
readsUsed=c(length(unlist(gread)))
cat("\n \n File contains",readsUsed,"reads that cover 128 bp windows of interest.\n",sep=" ")

if (max.reads<readsUsed) {cat("\n Only the first",max.reads,"reads for each region are used.")}

##Make smoothed averages for window

smoothWave=list()
### Allocate the memory required 
k = 1
ipd_nums = vector(length=length(gread),mode="integer")
for ( i in 1:length(gread)){
	idx = gread[[i]]
	if(length(idx) > 0){
	ipd = getIPD(h5,idx=idx)
	ipd_nums[i] = min(length(ipd),max.reads) 
	}
	else{
		ipd_nums[i] = 0
	}
}
#for ( i in 1:length(g0read)){
#	idx = g0read[[i]]
#	if(length(idx) > 0){
#	ipd = getIPD(h5,idx=idx)
#	ipd_nums[i+length(g1read)] = min(length(ipd),max.reads) 
#	}
#	else{
#		ipd_nums[i+length(g1read)] = 0
#	}
#}
no_reads = sum(ipd_nums) * 128
#print(no_reads)
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
##Remove inserts
	
	for ( i in 1:min(length(ipd),max.reads)){
		ins <- which(align[[i]][,2]=="-")
		
##If there are no inserts
  if(length(ins)==0) {
  	align[[i]]=align[[i]][,2]

if( total.time==TRUE) {
	instsTime=c(rep(0,length(ipd[[i]])));
	}

# If there are inserts
}else{
	if(total.time==TRUE){
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
		if (total.time==TRUE) instsTime <- instsTime[-ins];
		}
	if(reading_g0){
	  positz <- c(gR128[[k]][1]:gR128[[k]][2])	
	} else {
		positz <- c(gF128[[k-fwd_length]][2]:gF128[[k-fwd_length]][1])
	}
	align[[i]] <- align[[i]][match(positz,posit[[i]])]
	ipd[[i]] <- ipd[[i]][match(positz,posit[[i]])]
	pw[[i]] <- pw[[i]][match(positz,posit[[i]])]
	if (total.time==TRUE)	instsTime <- instsTime[match(positz,posit[[i]])];
	posit[[i]] <- positz
 	if(reading_g0){ 
	interp <- interp.0[posit[[i]]]
	}else{
	interp <- interp.1[posit[[i]]]
	}
	ipd[[i]][is.na(ipd[[i]])] <- 0
	pw[[i]][is.na(pw[[i]])] <- 0
	
	
	   if (total.time==TRUE){
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
	return(list(detailWave=detailWave,smoothWave=smoothWave,totalElements=totalElements,meanElementSize=meanElementSize,alignments=align,DNAPattern=DNAPattern))
})
	#########################	
			############################
		#########################	

#	if(length(g0read)!=0){
#	while (k <=length(g1read)+length(g0read)){
#		idx <- g0read[[k-length(g1read)]]
#			#print(i)
#		if(length(idx)>0){
#	posit <- getTemplatePosition(h5,idx=idx)
#	align <- getAlignments(h5,idx=idx)
#	ipd <- getIPD(h5,idx=idx)
#	pw <- getPulseWidth(h5,idx=idx)
#					
#	
#	for ( i in 1:min(length(ipd),max.reads)){
#		ins <- which(align[[i]][,2]=="-")#Find inserts
###If there are no inserts
#  if(length(ins)==0){
#  	 align[[i]]=align[[i]][,2]
#if (total.time==TRUE) instsTime=c(rep(0,length(ipd[[i]])));
#}
####If there are inserts
#	if(length(ins)>0){
#	if (total.time==TRUE){
#	instsTime=c(rep(0,length(ipd[[i]])))
#	tempITime=0
#	for ( j in 1:length(ins)){
#	  	pos=ins[j]
#		  nxtup=pos+1
#		  if(ins[j+1]>nxtup & j<length(ins)){
#		  	instsTime[nxtup] <- tempITime+ipd[[i]][pos]+pw[[i]][pos]
#			  tempITime=0
#		  }
#		  else{
#			  tempITime=tempITime+ipd[[i]][pos]+pw[[i]][pos]
#		  }
#	  }
#	}
###remove insert from all data		
#		
#		align[[i]] <- align[[i]][-ins,2]
#		ipd[[i]] <- ipd[[i]][-ins]
#		posit[[i]] <- posit[[i]][-ins]
#		pw[[i]] <- pw[[i]][-ins]
#		if (total.time==TRUE) instsTime <- instsTime[-ins];
#		}
#	
#		positz <- c(gR128[[ k - length(g1read)]][1]:gR128[[k - length(g1read) ]][2])
#	align[[i]] <- align[[i]][match(positz,posit[[i]])]
#	ipd[[i]] <- ipd[[i]][match(positz,posit[[i]])]
#	pw[[i]] <- pw[[i]][match(positz,posit[[i]])]
#	
#	if (total.time==TRUE)	instsTime <- instsTime[match(positz,posit[[i]])];
#
#	
#	posit[[i]] <- positz
#  
#	interp <- interp.0[posit[[i]]]
#
#	ipd[[i]][is.na(ipd[[i]])] <- 0
#	pw[[i]][is.na(pw[[i]])] <- 0
#		   
#		   if (total.time==TRUE){
#	   ipd.wst <- wd(rev(ipd[[i]]+pw[[i]]+instsTime),family="DaubExPhase",filter.number=filterNumber,type="station")
#   }
#   if (total.time==FALSE){
#   	ipd.wst <- wd(rev(ipd[[i]]),family="DaubExPhase",filter.number=filterNumber,type="station")
#	   }
#pat.wst <- wd(rev(interp),family="DaubExPhase",filter.number=filterNumber,type="station")
#	
#				st.index = 128 * (read_insert - 1) + 1
#				end.index = 128 * (read_insert)
#				#print(st.index)
#				#print(end.index)
#				#if(k == 1){
#				#	print(accessC(ipd.wst,level=7))
#				#}
#				#print(dim(smoothWave[[2]]))
#
#								for ( j in 1:8){
#							smoothWave[[j]][st.index:end.index,] = cbind(accessC(ipd.wst,level=j-1),accessC(pat.wst,level=j-1))
#								
#			}
#		if (plotDetailAvg==TRUE){
#	for (j in 1:3){
#			detailWave[[j]][st.index:end.index,] = matrix(accessD(ipd.wst,level=7-j),ncol=1)
#				}
#			}
#									
#	
#
#		read_insert = read_insert + 1	
#									
#		}
# }	
#	k=k+1
#	}
#}	
#if (plotSmoothAvg==T){
#			if(length(fileID)>0)pdf(paste(fileID,"smoothPlots.pdf",sep="_"))		
#		if(length(fileID)==0)pdf('smoothPlots.pdf')
#		
#		par(mfcol=c(3,1))	
#						
#mean0=c()
#qnt0=matrix(NA,nrow=4,ncol=128)
#
#						
#
#		for(i in 1:128){
#			idx=seq(i,length(smoothWave[[8]][,1]),by=128)
#				mean0[i]=mean(smoothWave[[8]][idx,1])
#qnt0[1,i]=quantile(smoothWave[[8]][idx,1],0.05)					
#qnt0[2,i]=quantile(smoothWave[[8]][idx,1],0.95)
#qnt0[3,i]=quantile(smoothWave[[8]][idx,1],0.1)
#qnt0[4,i]=quantile(smoothWave[[8]][idx,1],0.9)
#
#		}
#		
#		mainz=paste("IPD-wavelet for ",DNAPattern," (",nrow(smoothWave[[1]])/128,
#		" reads in ",totalElements," regions)",sep="")
#		poz=c()
#		if(totalElements==1){poz=positz[c(1,8*(1:16))]}else{poz=c(1,8*(1:16))}
#		if(shrink<1)(mainz=paste(shrink,"Shrunken mean IPD for ",DNAPattern," (",nrow(smoothWave[[8]])/128," reads/",totalElements,
#		" regions)",sep=""))
#		
#		plot(c(1:128),rev(mean0),ty="l",bty="n",ylab="IPD (s)",xlab="position in window",main=mainz,ylim=c(0,1.1*max(qnt0,mean0)),axes=F)
#		
#		axis(side=1,labels=poz,at=c(1,8*(1:16)))
#		axis(side=2)
#		
#		rect(shift.window+0.5,-100,shift.window+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
#		
#		segments(c(1:128),rev(qnt0[1,]),c(1:128),rev(qnt0[2,]),col=rgb(0,0,0,0.2))
#						segments(c(1:128),rev(qnt0[3,]),c(1:128),rev(qnt0[4,]),col=rgb(0,0,0,0.4))
#
#if (totalElements==1) text(1:128,1.05*(max(qnt0,mean0)),align[[1]],cex=0.5)
#	
#	if(totalElements!=1) {
#		text((shift.window-5),1.05*(max(qnt0,mean0)),"Read:",cex=0.8)
#		text((shift.window+1):(shift.window+floor(meanElementSize)),1.05*(max(qnt0,mean0)),align[[1]][(shift.window+1):(shift.window+floor(meanElementSize))],cex=0.5);}
#			
#
#				mean1=c()
#				qnt1=matrix(NA,nrow=4,ncol=128)
#				
#		for(i in 1:128){
#			idx=seq(i,length(smoothWave[[7]][,1]),by=128)
#				mean1[i]=mean(smoothWave[[7]][idx,1])
#qnt1[1,i]=quantile(smoothWave[[7]][idx,1],0.05)					
#qnt1[2,i]=quantile(smoothWave[[7]][idx,1],0.95)
#qnt1[3,i]=quantile(smoothWave[[7]][idx,1],0.1)
#qnt1[4,i]=quantile(smoothWave[[7]][idx,1],0.9)	
#							}
#		plot(c(1:128),rev(mean1),ty="l",bty="n",ylab=expression(E(s(2))),xlab="position in window",main="2 bp smoothing",ylim=c(0,1.1*max(qnt1,mean1)),axes=F)
#		
#		byz=c(8*(1:16))
#		
#		axis(side=1,labels=poz,at=c(1,byz))
#		axis(side=2)
#		
#		rect(shift.window+0.5,-100,shift.window+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
#			segments(c(1:128),rev(qnt1[1,]),c(1:128),rev(qnt1[2,]),col=rgb(0,0,0,0.2))
#						segments(c(1:128),rev(qnt1[3,]),c(1:128),rev(qnt1[4,]),col=rgb(0,0,0,0.4))
#
#				
#if (totalElements==1) text(1:128,1.05*(max(qnt1,mean1)),align[[1]],cex=0.5)
#
#	if(totalElements!=1) {
#		text((shift.window-5),1.05*(max(qnt1,mean1)),"Read:",cex=0.8)
#		text((shift.window+1):(shift.window+floor(meanElementSize)),1.05*(max(qnt1,mean1)),align[[1]][(shift.window+1):(shift.window+floor(meanElementSize))],cex=0.5);}
#			
#
#
#					mean2=c()
#					qnt2=matrix(NA,nrow=4,ncol=128)					
#		for(i in 1:128){
#			idx=seq(i,length(smoothWave[[6]][,1]),by=128)
#				mean2[i]=mean(smoothWave[[6]][idx,1])
#					
#qnt2[1,i]=quantile(smoothWave[[6]][idx,1],0.05)					
#qnt2[2,i]=quantile(smoothWave[[6]][idx,1],0.95)
#qnt2[3,i]=quantile(smoothWave[[6]][idx,1],0.1)
#qnt2[4,i]=quantile(smoothWave[[6]][idx,1],0.9)	
#
#
#		}
#		
#		plot(c(1:128),rev(mean2),ty="l",bty="n",ylab=expression(E(s(4))),xlab="position in window",main="4 bp smoothing",ylim=c(0,1.1*max(qnt2,mean2)),axes=F)
#		
#		axis(side=1,labels=poz,at=c(1,byz))
#		axis(side=2)
#		
#		rect(shift.window+0.5,-100,shift.window+meanElementSize+.5,100,col=rgb(0,0,1,0.1),border=NA)
#
#segments(c(1:128),rev(qnt2[1,]),c(1:128),rev(qnt2[2,]),col=rgb(0,0,0,0.2))
#						segments(c(1:128),rev(qnt2[3,]),c(1:128),rev(qnt2[4,]),col=rgb(0,0,0,0.4))
#				
#if (totalElements==1) text(1:128,1.05*(max(qnt2,mean2)),align[[1]],cex=0.5)	;
#	
#	if(totalElements!=1){
#				text((shift.window-5),1.05*(max(qnt2,mean2)),"Read:",cex=0.8)
#
#		 text((shift.window+1):(shift.window+floor(meanElementSize)),1.05*(max(qnt2,mean2)),align[[1]][(shift.window+1):(shift.window+floor(meanElementSize))],cex=0.5);}
#			
#			
#			
#		dev.off()	
#	}

################################################
################################################
################################################

#if (plotCorrelations==T & shrink==1){
#				if(length(fileID)>0) pdf(paste(fileID,"waveletCorrelations.pdf",sep="_"))
#				if(length(fileID)==0)pdf("waveletCorrelations.pdf")		
#		par(mfcol=c(1,2))	
#
#corz=list()
#pwz=matrix(NA,ncol=2,nrow=6)
#for(i in 1:6){
#corz[[i]]=cor.test(smoothWave[[i+1]][,1],smoothWave[[8]][,2])
##corz[[i]]=cor.test(smoothWave[[i+1]][,1],smoothWave[[i+1]][,2])
#pwz[i,1]=sum(smoothWave[[i+1]][,1]^2)
#pwz[i,2]=sum(smoothWave[[8]][,2]^2)
#}
#
#pwz[,1]=pwz[,1]/sum(pwz[,1])
#pwz[,2]=pwz[,2]/sum(pwz[,2])
#
#cor.val=rev(unlist(lapply(corz,"[",4)))
#cor.CI=rev(lapply(corz,"[",9))
#
#labels=c("2","4","8","16","32","64")
#
##ylimz=c(min(unlist(cor.CI)),
##max(unlist(cor.CI)))
#
#plot(cor.val,col="black",bty="n",axes=F,ylim=c(-1,1),xlab="scale (bp)",ylab="Pearson's r",main="Pattern-IPD Correlation",lwd=2,type="l")
#
#
#abline(h=c(-10:10)*0.1,col=rgb(0,0,0,0.1),lty=3)
#abline(h=0,col=rgb(0,0,0,0.2))
#
#for ( i in 1:6){
#lines(c(i,i),cor.CI[[i]]$conf.int[1:2],col="red",lty=1,lwd=2)
#}
#
#
#axis(1,at=1:6,labels=labels,las=1)
#axis(2,las=2,at=c(-10:10)*0.1,labels=F)
#axis(2,las=2,at=c(-1,-0.5,0,0.5,1),labels=T)
###$%%^&*
#plot(pwz[,1],ty='o',col='red',ylim=c(0,max(pwz)),main=paste("Power: IPD and", round(meanElementSize,1),"bp pattern"),xlab='scale',ylab='power',axes=F)
#points(pwz[,2],ty='o',col='blue')
#
#axis(1,at=1:6,labels=labels,las=1)
#axis(2,las=2,at=c(-10:10)*0.1,labels=F)
#axis(2,las=2,at=c(0:10)*0.1,labels=T)
#legend("topright",c("IPD","Pattern"),bty='n',col=c('red','blue'),lty=1)
#
#dev.off()
#}


#if (plotDetailAvg==T){
#			if(length(fileID)>0)pdf(paste(fileID,'detailPlots.pdf',sep="_"))
#			if(length(fileID)==0)pdf("detailPlots.pdf")		
#		
#		
#		par(mfcol=c(3,1))	
#						
#mean0=c()
#qnt0=matrix(NA,nrow=4,ncol=128)
#
#						
#
#		for(i in 1:128){
#			idx=seq(i,length(detailWave[[1]][,1]),by=128)
#				mean0[i]=mean(detailWave[[1]][idx,1])
#qnt0[1,i]=quantile(detailWave[[1]][idx,1],0.05)					
#qnt0[2,i]=quantile(detailWave[[1]][idx,1],0.95)
#qnt0[3,i]=quantile(detailWave[[1]][idx,1],0.1)
#qnt0[4,i]=quantile(detailWave[[1]][idx,1],0.9)
#
#		}
#		
#		mainz=paste("IPD-2bp detail wavelet for ",DNAPattern," (",nrow(detailWave[[1]])/128,
#		" reads in ",totalElements," regions)",sep="")
#		poz=c()
#		if(totalElements==1){poz=positz[c(1,8*(1:16))]}else{poz=c(0,8*(1:16))}
#		if(shrink<1)(mainz=paste(shrink,"Shrunken 2bp detail ",DNAPattern," (",nrow(detailWave[[1]])/128," reads/",totalElements,
#		" regions)",sep=""))
#		
#		plot(c(1:128),rev(mean0),ty="l",bty="n",ylab=expression(E(d(2))),xlab="position in window",main=mainz,ylim=c(min(qnt0,mean0),1.1*max(qnt0,mean0)),axes=F)
#		
#		axis(side=1,labels=poz,at=c(1,8*(1:16)))
#		axis(side=2)
#		
#		rect(shift.window+0.5,-100,shift.window+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
#
#		
#		segments(c(1:128),rev(qnt0[1,]),c(1:128),rev(qnt0[2,]),col=rgb(0,0,0,0.2))
#						segments(c(1:128),rev(qnt0[3,]),c(1:128),rev(qnt0[4,]),col=rgb(0,0,0,0.4))
#
#if (totalElements==1) text(1:128-0.5,1.05*(max(qnt0,mean0)),align[[1]],cex=0.5)
#
#if(totalElements!=1){
#			text((shift.window-5),1.05*(max(qnt0,mean0)),"Read:",cex=0.8)
#
#	 text((shift.window+1):(shift.window+floor(meanElementSize)),1.05*(max(qnt0,mean0)),align[[1]][(shift.window+1):(shift.window+floor(meanElementSize))],cex=0.5);}
#			
#			
#
#				mean1=c()
#				qnt1=matrix(NA,nrow=4,ncol=128)
#				
#		for(i in 1:128){
#			idx=seq(i,length(detailWave[[2]][,1]),by=128)
#				mean1[i]=mean(detailWave[[2]][idx,1])
#qnt1[1,i]=quantile(detailWave[[2]][idx,1],0.05)					
#qnt1[2,i]=quantile(detailWave[[2]][idx,1],0.95)
#qnt1[3,i]=quantile(detailWave[[2]][idx,1],0.1)
#qnt1[4,i]=quantile(detailWave[[2]][idx,1],0.9)	
#							}
#		plot(c(1:128),rev(mean1),ty="l",bty="n",ylab=expression(E(d(4))),xlab="position in window",main="4 bp detail",ylim=c(min(qnt1,mean1),1.1*max(qnt1,mean1)),axes=F)
#		
#		byz=c(8*(1:16))
#		
#		axis(side=1,labels=poz,at=c(1,byz))
#		axis(side=2)
#		
#		rect(shift.window+0.5,-100,shift.window+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
#			segments(c(1:128),rev(qnt1[1,]),c(1:128),rev(qnt1[2,]),col=rgb(0,0,0,0.2))
#						segments(c(1:128),rev(qnt1[3,]),c(1:128),rev(qnt1[4,]),col=rgb(0,0,0,0.4))
#
#				
#if (totalElements==1) text(1:128-0.5,1.05*(max(qnt1,mean1)),align[[1]],cex=0.5)
#
#if(totalElements!=1) {
#			text((shift.window-5),1.05*(max(qnt1,mean1)),"Read:",cex=0.8)
#
#	text((shift.window+1):(shift.window+floor(meanElementSize)),1.05*(max(qnt1,mean1)),align[[1]][(shift.window+1):(shift.window+floor(meanElementSize))],cex=0.5);}
#			
#			
#
#					mean2=c()
#					qnt2=matrix(NA,nrow=4,ncol=128)					
#		for(i in 1:128){
#			idx=seq(i,length(detailWave[[3]][,1]),by=128)
#				mean2[i]=mean(detailWave[[3]][idx,1])
#					
#qnt2[1,i]=quantile(detailWave[[3]][idx,1],0.05)					
#qnt2[2,i]=quantile(detailWave[[3]][idx,1],0.95)
#qnt2[3,i]=quantile(detailWave[[3]][idx,1],0.1)
#qnt2[4,i]=quantile(detailWave[[3]][idx,1],0.9)	
#
#
#		}
#		
#		plot(c(1:128),rev(mean2),ty="l",bty="n",ylab=expression(E(d(8))),xlab="position in window",main="8 bp detail",ylim=c(min(qnt2,mean2),1.1*max(qnt2,mean2)),axes=F)
#		
#		axis(side=1,labels=poz,at=c(1,byz))
#		axis(side=2)
#		
#		rect(shift.window+0.5,-100,shift.window+meanElementSize+0.5,100,col=rgb(0,0,1,0.1),border=NA)
#
#segments(c(1:128),rev(qnt2[1,]),c(1:128),rev(qnt2[2,]),col=rgb(0,0,0,0.2))
#						segments(c(1:128),rev(qnt2[3,]),c(1:128),rev(qnt2[4,]),col=rgb(0,0,0,0.4))
#				
#if (totalElements==1) text(1:128-0.5,1.05*(max(qnt2,mean2)),align[[1]],cex=0.5)	
#		
#		if(totalElements!=1) {
#					text((shift.window-5),1.05*(max(qnt2,mean2)),"Read:",cex=0.8)
#
#			text((shift.window+1):(shift.window+floor(meanElementSize)),1.05*(max(qnt2,mean2)),align[[1]][(shift.window+1):(shift.window+floor(meanElementSize))],cex=0.5);}
#			
#			
#			
#		dev.off()	
#	}
#
#
#}

#End script;
