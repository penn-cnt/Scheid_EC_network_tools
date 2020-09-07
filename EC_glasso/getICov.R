library(bootnet)
library(R.matlab)

testmat<-readMat('tempSet.mat')
data=testmat$data
gamma=testmat$gamma
N=dim(data)[1]
networks<-array(NA,c(N,N,testmat$TT))
lambdas<-array(NA, c(testmat$TT))
icovRho<-array(NA,c(N,N,testmat$TT))

for (i in 1:testmat$TT){
  cat(i, "\n")
  window<- data[seq(1,N),seq((i-1)*testmat$Lwin+1,i*testmat$Lwin)]
  results <- estimateNetwork(
    t(window),
    default = "EBICglasso",
    tuning = gamma,
    corMethod = "npn",
    lambda.min.ratio = 0.01,
    penalize.diagonal = TRUE)
  
  l=1
  
  # Get index for the optimal tuning matrix
  idx = which(results$results$optwi[l,l] == results$results$results$wi[1,1,])
  #idx2= which(results$results$lambda >= 0.01)
  
  networks[,,i] = results$graph
  lambdas[i]=results$results$lambda[idx]
  icovRho[,,i]=results$results$optwi
  
}

writeMat('glasso.mat',nets=networks,lambdas=lambdas, icovRho=icovRho)

# #View image
# 
# penDiag=TRUE
#  #localresults<-glasso(cov(t(window)), 0.01, penalize.diagonal=penDiag)
#  localresults<-glasso(cc$cc), 0.01, penalize.diagonal=penDiag)
#  pathresults<-glassopath(cov(t(window)), 0.01, penalize.diagonal=penDiag, trace=0)
#  writeMat('Rlasso.mat', Rlasso=localresults, penDiag=penDiag, glassoPath=pathresults$wi)
