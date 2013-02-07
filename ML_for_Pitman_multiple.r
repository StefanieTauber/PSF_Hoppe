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

## param: [type = "vector"]: 
## depicts the parameter vector c(theta,sigma)

## freq: [type = "list"]:
## each slot of the list contains the number of reads starting at each position for one gene (G-level case). Positions where 0 reads start can be neglected. The order of the slots does not matter. The length of the list corresponds to the number of genes you want to integrate in the analysis

## k: [type = "numeric vector"]:
## contains the number of positions where at least 1 read starts for each gene (G-level case) 

## n: [type = "numeric vector"]: 
## contains the sum of reads mapping to the respective gene (G-level case)


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
 
 
############################################################################################ 
### PSF for finite sampling universe, optimize parameters for multiple frequency vectors ###
############################################################################################


## finite sampling universe: theta = m*|sigma| and sigma is smaller zero
## in the G-level case m = length of the gene (not more than "length of gene" positions can possibly get selected as SP
## Log of Pitman is easier to compute
## freq, k, n: like above

## lot: [type = "numeric vector"]
## contains the length of the corresponding genes


LogPitman_multiple_finite <- function(theta,freq,k,n,lot) {

  tmp = 0
  
  for (i in 1:length(k)){
  
  tmp = tmp + sum(log(theta + c(1:(k[i]-1))*(-theta/lot[i])))- LogAscFact(theta+1,n[i]-1)+sum(mapply(LogAscFact,rep(1-(-theta/lot[i]),k[i]),freq[[i]]-1))}
  
  return(tmp)}
  
  
## example use: 
# optimize(LogPitman_multiple_finite, interval = c(1,10000), freq, k, n,lot)$maximum
  

