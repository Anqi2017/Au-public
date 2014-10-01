####### entropy_reduce.R ##########
# Entropy reduce reduces a set of genes representing a group of experiments by removing
# genes that contribute the least to the important eigenvectors.  It reduces the number
# of genes one by one until the entropy of the system is reduced by some fraction.  The
# number of important eigenvectors to be considered is determined by a fraction of
# variance explained by the eigenvectors
#
# Pre: A filename for an input table of genes by row for experiments by columns, with unique labels for each
#      A filename for an output table
#      Hardcoded:  input tolerance (default = 0.1) for a number considered close to zero for rankMatrix
#                  input threshold (for entropy) (default 0.9) amount the entropy should be reduced before exiting
#                        this is a little funky but on datasets with mostly very informative genes, removing the least informative genes increases entropy.
#                        therefor we also exit if you move the same distance in the other way. (so if 0.9 is input, then 0.9 1.1 will be the bounds).
#                  input variance threshold (default 0.95) consider the contributions to eigenvectors out to a number that explains this percent of the variance. 0 would be the first eigenvector only, 1 would be all eigenvectors.
# Post: A table formatted exactly like the input minus the removed rows
# Modifies: STDOUT shows progress
#
library('methods')
library('lattice')
library('Matrix')
args <- commandArgs(trailingOnly = TRUE)
inputfile = args[1]
outputfile = args[2]
table = read.table(inputfile,header=TRUE,row.names=1, sep="\t",check.names=FALSE)
inputtol = 0.1
inputthresh = 0.90
inputvarthresh = 0.95
rowi = row.names(table)
si = svd(table)
#calculate the relative contribution of each eigenvector
pi = si$d/sum(si$d)
vi = (si$d*si$d)/sum(si$d*si$d)
#calculate the entropy (Alter et. al. PNAS 2000)
ranki = rankMatrix(table,tol=inputtol)[1]
ei = -1/log(ranki)*sum((vi*log(vi))[1:ranki])
ethresh = ei*inputthresh
ethresh2 = ei*((1-inputthresh)+1)
p = pi # start our looping with the initial weights of the eigenvectors
v = vi
s = si # start our looping the the intial svd
e = ei # start looping with the intial entropy
runtable = table # start our looping with the intial table
while(e > ethresh && e < ethresh2) {
  savetable = runtable
  #abs here could be contentious as it counts either positive or negative as important.
  # also not sure whether to use v or p the variance proportion or the standard deviation proportion.
  contrib = abs(s$u) %*% diag(v)
  var = 0
  #cnt will be the number of eigenvectors to use based on the
  #    percent of variance they explain
  cnt = 0
  while(var < inputvarthresh && cnt < length(v)) {
    cnt = cnt + 1
    var = sum(v[1:cnt])
  }
  rowval = rowSums(contrib[,1:cnt])
  worstrowindex = which.min(rowval)
  print(cbind(dim(runtable)[1],row.names(runtable)[worstrowindex],rowval[worstrowindex],e,ethresh,cnt))
  if(worstrowindex == 1) {
    runtable = runtable[2:dim(runtable)[1],]
  } else if(worstrowindex == dim(s$u)[1]) {
    runtable = runtable[1:(worstrowindex-1),]
  } else {
    runtable = rbind(runtable[1:(worstrowindex-1),],runtable[(worstrowindex+1):dim(runtable)[1],])
  } 
  s = svd(runtable)
  p = s$d/sum(s$d) #proportion of standard deviation
  v = (s$d*s$d)/sum(s$d*s$d) #proportion of variance
  #e = -1/log(rankMatrix(runtable,tol=inputtol)[1])*sum(v*log(v))
  rnk = rankMatrix(runtable,tol=inputtol)[1]
  e = -1/log(rnk)*sum((v*log(v))[1:rnk])
  #print (e)
}
write.table(savetable,file=outputfile,sep="\t",quote=FALSE)
