#Sebas Ramos-Onsins 18/Feb/2014

args <- commandArgs(trailingOnly=TRUE) #collect data from arguments
#first argument: nchromosomes
#second argument: file with the total length of each chromosome, on rows.
#third argument: file with the name of marker, number of chromosome, cM and Physical position (in bp)
#fourth and rest arguments: size or sizes (in bp) of the windows.

#function derivative
derivative.triangle <- function(f1,f2,p1,p2) {
	slp <- (f2-f1)/(p2-p1)
	return(slp)
}

#nchr <- 12
#len.chr <- array(c(35383099,26193771,29387469,33123231,28337775,35939859,26773857,32513408,24107567,25362316, 31442130,26400393),dim=c(nchr)) #the length of each chromosome
#RecMap <- read.table(file="marker_all_LG.txt",header=TRUE)
#window.sizes <- c(5e4,1e5,5e5,1e6,5e6)

nchr <- as.numeric(args[1])
len.chr <- read.table(file=args[2],header=TRUE)
len.chr <- as.matrix(len.chr)
RecMap <- read.table(file=args[3],header=TRUE)
window.sizes <- c(as.numeric(args[4:length(args)]))

#show(nchr)
show(len.chr)
#show(RecMap)
#show(window.sizes)
show(length(len.chr))

len.chr.cum <- len.chr
for(i in 2:length(len.chr)) {
	len.chr.cum[i] <- len.chr.cum[i] + len.chr.cum[i-1]
}
len.chr.cum <- c(0,len.chr.cum)


RecMap <- RecMap[,c(2,3,4)]
head(RecMap)

#include the first bp position
initpositions <- data.frame(cbind(LG=c(1:nchr),cM=rep(0,nchr),ChrPos=rep(1,nchr)))
RecMap  <-  rbind(RecMap,as.data.frame(initpositions))
RecMap <- RecMap[order(RecMap[,1],RecMap[3]),] #sort using the first and the third columns
head(RecMap)

#erasing those positions with cM smaller than previous: (...)
RecMap2 <- RecMap
i <- 2
while(i<length(RecMap2[,1])) {
	if((RecMap2[i,1] == RecMap2[i-1,1]) && (RecMap2[i,2] < RecMap2[i-1,2])) {
		RecMap2 <- RecMap2[-i,]
	} else {
		i <- i + 1
	}
}
dim(RecMap)
dim(RecMap2)

for(window in window.sizes) {
	rec.win <- data.frame(cbind(chr=0,start=0,end=0,len=0,recw=0))
	nwin <- 1
	
	pdf(sprintf("RecMap_%d.pdf",window))
	par(mar=c(5,4,4,5))
	for(chr in 1:nchr) {
		#do the spline:
		ww <- array(1,dim=c(floor(len.chr[chr]/window)+1))
		for(y in 2:length(ww)) {ww[y] <- ww[y-1] + window}
		last.position <- RecMap2[RecMap2[,1]==chr,3][length(RecMap2[RecMap2[,1]==chr,3])]
		ww <- ww[ww < last.position]
		ww <- c(ww,last.position)
		cum.rec <- spline(RecMap2[RecMap2[,1]==chr,3],RecMap2[RecMap2[,1]==chr,2],method="hyman",xout=ww)

		#plot the cummulative curves
		plot(RecMap2[RecMap2[,1]==chr,3],RecMap2[RecMap2[,1]==chr,2],main=sprintf("Cummulative Recombination Map\n Chr%d Windows of %0.f bp",chr,window),xlab=sprintf("Physical Distance (bp)"),ylab=sprintf("Genetical Distance (cM)"))
		lines(cum.rec,col="blue")
		
		#calculate the recombination per window
		nwin.start <- nwin
		chrd  <- sprintf("chr%02d",chr)
		for(w in 1:(length(cum.rec$x)-1)) {
			recw <- derivative.triangle(cum.rec$y[w],cum.rec$y[w+1],cum.rec$x[w],cum.rec$x[w+1])
			if(recw < 0) {recw <- 0}
			rec.win[nwin,] <- c(chr,as.numeric(sprintf("%.0f",cum.rec$x[w])),as.numeric(sprintf("%.0f",cum.rec$x[w+1]-1)),as.numeric(sprintf("%.0f",cum.rec$x[w+1]-cum.rec$x[w])),recw)
			rec.win <- rbind(rec.win,cbind(chr=0,start=0,end=0,len=0,recw=0))
			nwin <- nwin + 1
		}
		par(new=T)
		plot(cum.rec$x[-1]-1/2*window,rec.win$recw[c(nwin.start:(nwin-1))],ylim=c(0,2*max(rec.win$recw[c(nwin.start:nwin)])),axes=F,lty=1,col="red",ylab="",xlab="",type="l")
		axis(side=4)
		mtext(4,text="Recombination Rate",line=3)
	}
	dev.off() #close the plot file
	
	#print and draw recombination per position
	rec.win <- rec.win[-length(rec.win[,1]),]
	rec.win

	start_chr <- c(which(rec.win[,2]==1),length(rec.win[,2]))
	positions <- array(rec.win[,2],dim=c(length(rec.win[,1])))
	for(w in 1:length(rec.win[,2])) {
		positions[w] <- positions[w] + len.chr.cum[rec.win[w,1]]
	}
	start_chr <- positions[start_chr]
	write.table(cbind(rec.win,positions),quote=FALSE,sep="\t",row.name=FALSE,file=sprintf("RecMap_windows_%d.txt", window))

	#include one more row to start each plot by NA recw
	initpositions <- data.frame(cbind(chr=c(1:nchr),start=rep(0,nchr),end=rep(0,nchr),len=rep(0,nchr),recw=rep(NA,nchr)))
	rec.win2  <-  rbind(rec.win,as.data.frame(initpositions))
	rec.win2 <- rec.win2[order(rec.win2[,1],rec.win2[2]),] #sort using the first and the second columns
	start_chr <- c(which(rec.win2[,2]==1),length(rec.win2[,2]))
	positions <- array(rec.win2[,2],dim=c(length(rec.win2[,1])))
	for(w in 1:length(rec.win2[,2])) {
		if(positions[w] > 1) positions[w] <- positions[w] - 1/2*window
		positions[w] <- positions[w] + len.chr.cum[rec.win2[w,1]]
	}
	start_chr <- positions[start_chr]
	
	pdf(sprintf("RecMap_all_chromosomes_%d.pdf",window),width=60,height=10)
	par(mfrow=c(1,1))
	plot(positions,rec.win2$recw,type="l",lty=1,pch=20,cex=0.5,col="black",main=sprintf("Recombination map using window size of %d",window),xlab="Positions",ylab=colnames(rec.win)[5])
	lines(abline(v=start_chr,col="grey",lty=1))
	dev.off()	
}
