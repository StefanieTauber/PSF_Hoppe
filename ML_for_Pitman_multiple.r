####################################
### ascending factorial function ###
####################################

LogAscFact <- function(a,n) {
  if (n > 0) {
    return(sum(log(a+c(0:(n-1)))))
  } else {
    return(0)
  }
}

############################################################################################## 
### PSF for infinite sampling universe, optimize parameters for multiple frequency vectors ###
##############################################################################################

## typical use-case: if you want to estimate one common theta and sigma for distribution of reads within a gene over all genes

## param: [type = "vector"]: 
## depicts the parameter vector c(theta,sigma)

## freq: [type = "list"]:
## each slot of the list contains the number of reads starting at each position for one gene. Positions where 0 reads start can be neglected. The order of the slots does not matter. The length of the list corresponds to the number of genes you want to integrate in the analysis

## k: [type = "numeric vector"]:
## contains the number of positions where at least 1 read starts for each gene 

## n: [type = "numeric vector"]: 
## contains the sum of reads mapping to the respective gene


## infinite sampling universe: sigma within [0,1) and theta > - sigma
## Log of Pitman is easier to compute

LogPitman_multiple_infinite <- function(param,freq,k,n) {
     
     if (param[2] < 0 | param[2] >= 1 | param[1] <= -param[2] ) {
      return(NA)}
      
     else {
     tmpPitman = 0
    
    for (i in 1:length(k)){
    
      if (k[i] == 2 & param[1] < 0) {return(NA)}
      else {
    
      tmpPitman <- tmpPitman + sum(log(param[1] + c(1:(k[i]-1))*param[2]))-LogAscFact(param[1]+1,n[i] -1)+sum(mapply(LogAscFact,rep(1-param[2],k[i]),freq[[i]]-1))}
	i = i+1
    }}
    return(tmpPitman)
}

Gradient_LogPitman_multiple_infinite <- function(param,n,k,sum_n) {

    if (param[2] < 0 | param[2] >= 1 | param[1] <= -param[2]) {
    return(NA)}

    else {
    tmp = c(0,0)
      
       for (i in 1:length(k)){
       
        tmp[1] <- tmp[1] + sum(1/(param[1] + c(1:(k[i]-1))*param[2])) - sum(1/(param[1]+c(1:(n[i]-1))))
    
    tmp[2] <- tmp[2] + sum(c(1:(k[i]-1))/(param[1] + c(1:(k[i]-1))*param[2]))
    
    for(j in which(freq[[i]]>1)) {
      tmp[2] <- tmp[2] + sum(1/(param[2]-c(1:(freq[[i]][j]-1))))
    }}
    i = i+1}
    return(tmp)
  }

## example use: 
#  optim(par=c(0.01,0.01), fn = LogPitman_multiple_infinite, gr = Gradient_LogPitman_multiple_infinite, freq = freq, k = k, n = n, control = list(fnscale=-1)
 
 


