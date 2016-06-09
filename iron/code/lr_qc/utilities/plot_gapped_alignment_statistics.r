#!/opt/R/3.2.1/bin/Rscript

args=commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Must supply input and output files\n",call.=FALSE)
}
filex = substr(args[2],nchar(args[2])-2,nchar(args[2]))

if(filex=="pdf") {
  pdf(args[2])
} else if(filex=="png") {
  png(args[2],bg="#00000000")
} else {
  stop("Unsupported type of output file\n",call.=FALSE)
}  

infilex = substr(args[1],nchar(args[1])-1,nchar(args[1]))
if(infilex=="gz") {
  d<-read.table(args[1])
} else {
  d<-read.table(gzfile(args[1],open="r"))
}

fname = '/Shared/Au/jason/Code/NEWFUZZ/tempall/data/lengths.txt.gz'
d<-read.table(gzfile(fname,open="r"))
par(bg="#FFFFFF")
par(mar=c(4,4,1.5,0.5),oma=c(0.5,1,0.5,1))
mat=rbind(c(1,2,3),c(4,5,6))
layout(mat,c(1,5),c(6,6,3))

##### Pie chart for split alignments ##########

recwid = 1

multi = length(d[d[,2]!="original" & d[,2]!="unaligned",1])
single = length(d[d[,2]=="original",1])
unaligned = length(d[d[,2]=="unaligned",1])
tot = multi+single+unaligned
multi_perc = paste(round(100*multi/tot,1),"%",sep='')
single_perc = paste(round(100*single/tot,1),"%",sep='')
unaligned_perc = paste(round(100*unaligned/tot,1),"%",sep='')
acols = c("#FFFFFF","#777777","#FF0000")
alabels=c(paste("Unaligned",unaligned_perc),paste("Single best",single_perc),paste("Split best",multi_perc))
plot(1,type="n",xlim=c(0,500),ylim=c(0,tot+100),xaxt='n',ylab="Counts any length read",bty="n",xlab="",cex.axis=1.2,cex.lab=1.2)
rect(0,0,500,single,col="#777777",lwd=recwid)
rect(0,single,500,single+multi,col="#FF0000",lwd=recwid)
rect(0,single+multi,500,single+multi+unaligned,col="#FFFFFF",lwd=recwid)
mtext("All reads",side=1,at=-100,adj=0)


###### Read length-wise bar plot for split alignmetns #####
longest = 4000

#find the biggest bin
biggest = 0
for(i in seq(0,longest-500,500)) {
  currsize = length(d[d[,4]>i & d[,4]<=i+500,1])
  if(currsize > biggest) { biggest = currsize }
}
currsize = length(d[d[,4]>longest,1])
if(currsize > biggest) { biggest = currsize }


plot(1,type="n",xlim=c(0,longest),ylim=c(0,biggest+100),ylab="Count reads per length",bty="n",xlab="Read length (bp)",xaxt='n',cex.axis=1.2,cex.lab=1.2)
axis(side=1,at=seq(0,longest,500),cex.axis=1.2)
for(i in seq(0,longest-500,500)) {
  single = length(d[d[,2]=="original" & d[,5]>i & d[,5]<=i+500,1])
  rect(i,0,i+500,single,col="#777777",lwd=recwid)
  multi = length(d[d[,2]!="original" & d[,2]!="unaligned"  & d[,5]>i & d[,5]<=i+500,1])
  rect(i,single,i+500,single+multi,col="#FF0000",lwd=recwid)
  unalign = length(d[d[,2]=="unaligned" & d[,5]>i & d[,5]<=i+500,1])
  rect(i,single+multi,i+500,single+multi+unalign,col="#FFFFFF",lwd=recwid)
}

###last case###
single = length(d[d[,2]=="original" & d[,5]>longest,1])
multi = length(d[d[,2]!="unaligned" & d[,2]!="original" & d[,5]>longest,1])
unalign = length(d[d[,2]=="unaligned" & d[,5]>longest,1])
tot = single+multi+unalign
plot(1,type="n",xlim=c(0,500),ylim=c(0,tot+100),ylab="Count longest reads",bty="n",xlab="",xaxt='n',cex.axis=1.2,cex.lab=1.2)
rect(0,0,500,single,col="#777777",lwd=recwid)
rect(0,single,500,single+multi,col="#FF0000",lwd=recwid)
rect(0,single+multi,500,single+multi+unalign,col="#FFFFFF",lwd=recwid)
mtext(paste(">",longest),side=1,at=-100,adj=0)

plot(1,type='n',xlim=c(0,500),ylim=c(0,1),bty="n",xaxt='n',ylab="Fraction aligned any length read",xlab="",cex.axis=1.2,cex.lab=1.2)
dat1 = d[d[,2]!="unaligned",4]
dat2 = d[d[,2]!="unaligned",5]
dcom = dat1/dat2
boxplot(dcom,add=TRUE,at=250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE)
mtext("All reads",side=1,at=-100,adj=0)


####### Plot the distribution of reads ###########
plot(1,type='n',xlim=c(0,longest),ylim=c(0,1),bty="n",xaxt='n',ylab="Fraction of reads aligned",xlab="Read length(bp)",cex.axis=1.2,cex.lab=1.2)
axis(side=1,at=seq(0,longest,500),cex.axis=1.2)
points(d[d[,2]=="original",5],d[d[,2]=="original",4]/d[d[,2]=="original",5],col="#77777730",pch='.')
points(d[d[,2]!="original" & d[,2]!="unaligned",4],d[d[,2]!="original" & d[,2]!="unaligned",4]/d[d[,2]!="original" & d[,2]!="unaligned",5],col="#FF000030",pch='.')
for(i in seq(0,longest-500,500)) {
  dat1 = d[d[,2]!="unaligned" & d[,5]>i & d[,5]<=i+500,4]
  dat2 = d[d[,2]!="unaligned" & d[,5]>i & d[,5]<=i+500,5]
  dcom = dat1/dat2
  boxplot(dcom,add=TRUE,at=i+250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE)
}
dat1 = d[d[,2]!="unaligned" & d[,5]>longest,4]
dat2 = d[d[,2]!="unaligned" & d[,5]>longest,5]
dcom = dat1/dat2
plot(1,type='n',xlim=c(0,500),ylim=c(0,1),bty="n",xaxt='n',ylab="Fraction of longest reads aligned",xlab="",cex.axis=1.2,cex.lab=1.2)
boxplot(dcom,add=TRUE,at=250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE,col="#00000000")
mtext(paste(">",longest),side=1,at=-100,adj=0)

# ### add legend information to bottom
# plot(1,type='n',xlim=c(0,1),ylim=c(0,1),bty="n",xaxt='n',yaxt='n',ylab="",xlab="")
# legend(0.7,1,c(paste("Unaligned"),paste("Gapped Alignment"),paste("Single Alignment")),fill=c("#FFFFFF","#FF0000","#777777"),bty='n',cex=1.2)

dev.off()
