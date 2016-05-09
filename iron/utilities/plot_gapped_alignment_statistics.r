#!/opt/R/3.2.1/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)<2) {
  stop("Must supply input and output files\n",call.=FALSE)
}
pdf(args[2])
d<-read.table(args[1])

mat=cbind(c(1,3),c(2,4))
layout(mat,c(1.3,4),c(4,4))

##### Pie chart for split alignments ##########

multi = length(d[d[,1]>1,1])
single = length(d[d[,1]==1,1])
unaligned = length(d[d[,1]==0,1])
tot = multi+single+unaligned
multi_perc = paste(round(100*multi/tot,1),"%",sep='')
single_perc = paste(round(100*single/tot,1),"%",sep='')
unaligned_perc = paste(round(100*unaligned/tot,1),"%",sep='')
acols = c("#FFFFFF","#777777","#FF0000")
alabels=c(paste("Unaligned",unaligned_perc),paste("Single best",single_perc),paste("Split best",multi_perc))
plot(1,type="n",xlim=c(0,500),ylim=c(0,tot+100),xaxt='n',ylab="Counts any length read",bty="n",xlab="")
rect(0,0,500,single,col="#777777")
rect(0,single,500,single+multi,col="#FF0000")
rect(0,single+multi,500,single+multi+unaligned,col="#FFFFFF")


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


plot(1,type="n",xlim=c(0,longest+500+400),ylim=c(0,biggest+biggest*0.25),ylab="Count reads per length",bty="n",xlab="Read length (bp)",xaxt='n')
axis(side=1,at=seq(0,longest,500))
mtext(paste(">",longest),side=1,at=longest+400,adj=0)
for(i in seq(0,longest-500,500)) {
  single = length(d[d[,1]==1 & d[,4]>i & d[,4]<=i+500,1])
  rect(i,0,i+500,single,col="#777777")
  multi = length(d[d[,1]>1 & d[,4]>i & d[,4]<=i+500,1])
  rect(i,single,i+500,single+multi,col="#FF0000")
  unalign = length(d[d[,1]==0 & d[,4]>i & d[,4]<=i+500,1])
  rect(i,single+multi,i+500,single+multi+unalign,col="#FFFFFF")
}
#last case
single = length(d[d[,1]==1 & d[,4]>longest,1])
rect(longest+400,0,longest+400+500,single,col="#777777")
multi = length(d[d[,1]>1 & d[,4]>longest,1])
rect(longest+400,single,longest+400+500,single+multi,col="#FF0000")
unalign = length(d[d[,1]==0 & d[,4]>longest,1])
rect(longest+400,single+multi,longest+400+500,single+multi+unalign,col="#FFFFFF")
legend(1500,biggest+biggest*0.25,c(paste("Unaligned",unaligned_perc),paste("Gapped Alignment",multi_perc),paste("Single Alignment",single_perc)),fill=c("#FFFFFF","#FF0000","#777777"))

plot(1,type='n',xlim=c(0,500),ylim=c(0,1),bty="n",xaxt='n',ylab="Fraction aligned any length read",xlab="")
dat1 = d[d[,1]>=1,3]
dat2 = d[d[,1]>=1,4]
dcom = dat1/dat2
boxplot(dcom,add=TRUE,at=250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE)


####### Plot the distribution of reads ###########
plot(1,type='n',xlim=c(0,longest+400+500),ylim=c(0,1),bty="n",xaxt='n',ylab="Fraction of read aligned",xlab="Read length(bp)")
axis(side=1,at=seq(0,longest,500))
mtext(paste(">",longest),side=1,at=longest+400,adj=0)
points(d[d[,1]==1,4],d[d[,1]==1,3]/d[d[,1]==1,4],col="#77777705",pch='*')
points(d[d[,1]>1,4],d[d[,1]>1,3]/d[d[,1]>1,4],col="#FF000005",pch='*')
for(i in seq(0,longest-500,500)) {
  dat1 = d[d[,1]>=1 & d[,4]>i & d[,4]<=i+500,3]
  dat2 = d[d[,1]>=1 & d[,4]>i & d[,4]<=i+500,4]
  dcom = dat1/dat2
  boxplot(dcom,add=TRUE,at=i+250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE,outline=FALSE)
}
dat1 = d[d[,1]>=1 & d[,4]>longest,3]
dat2 = d[d[,1]>=1 & d[,4]>longest,4]
dcom = dat1/dat2
boxplot(dcom,add=TRUE,at=longest+400+250,boxwex=900,border="#000000",frame=FALSE,axes=FALSE)
dev.off()
