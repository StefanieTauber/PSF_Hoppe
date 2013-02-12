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


############################################################################################### 
### PSF for infinite sampling universe, optimize parameters for one single frequency vector ###
###############################################################################################

## param: [type = "vector"]: 
## depicts the parameter vector c(theta,sigma)

## freq: [type = "numeric vector"]:
## give the read counts per gene or the number of reads starting at each position. Genes with 0 reads or positions where 0 reads start can be neglected. The order of freq does not matter

## k: [type = "numeric"]:
## number of genes with at least 1 count or number of positions where at least 1 read starts 

## n: [type = "numeric"]: 
## sum of counts or sum of reads mapping to the gene 

## infinite sampling universe: sigma within [0,1) and theta > - sigma
## Log of Pitman is easier to compute

LogPitman_single_infinite <- function(param,freq,k,n) {
  if (param[2] < 0 | param[2] >= 1 | param[1] <= -param[2]) {
    return(NA)
  } else {
    return(sum(log(param[1] + c(1:(k-1))*param[2]))-LogAscFact(param[1]+1,n -1)+sum(mapply(LogAscFact,rep(1-param[2],k),freq-1)))
  }
}

Gradient_LogPitman_single_infinite <- function(param,freq,k,n) {
  if (param[2] < 0 | param[2] >= 1 | param[1] <= -param[2]) {
    return(NA)
  } else {
    d <- c(0,0)
    d[1] <- sum(1/(param[1] + c(1:(k-1))*param[2])) - sum(1/(param[1]+c(1:(n-1))))
    d[2] <- sum(c(1:(k-1))/(param[1] + c(1:(k-1))*param[2]))
    for(j in which(freq>1)) {
      d[2] <- d[2] + sum(1/(param[2]-c(1:(freq[j]-1))))
    }
    return(d)
  }
}


## example use: 
#  optim(par=c(0.01,0.01), fn = LogPitman_single_infinite, gr = Gradient_LogPitman_single_infinite, freq = freq, k = k, n = n, control = list(fnscale=-1))


 
############################################################################################# 
### PSF for finite sampling universe, optimize parameters for one single frequency vector ###
#############################################################################################


## finite sampling universe: theta = m*|sigma| and sigma is smaller zero
## m might be either the length of the gene or the maximum number of transcripts that may be expressed in the respective sequencing sample

## freq, k, n: like above

## lot: [type = "numeric"] size of transcript universe and, respectivley, length of respective transcript or gene



LogPitman_single_finite <- function(theta,freq,k,n,lot) {
      return(sum(log(theta + c(1:(k-1))*(-theta/lot)))-
      LogAscFact(theta+1,n-1)+sum(mapply(LogAscFact,rep(1-(-theta/lot),k),freq-1)))}
      


## example use: 
# optimize(LogPitman_single_finite, interval = c(1,10000), freq, k, n,lot,maximum = T)$maximum


###################################################################################################
###                                                                                             ###          
### PSF for finite sampling universe, optimize parameters for one single frequency vector       ###
### Here we are not optimizing theta and sigma but theta and m ( = extent of sampling universe) ###
###  (theta = m|sigma| holds)                                                                   ###    
###################################################################################################

## param: [type = "vector"]: 
## depicts the parameter vector c(theta,m)

## freq: [type = "numeric vector"]:
## give the read counts per gene or the number of reads starting at each position. Genes with 0 reads or positions where 0 reads start can be neglected. The order of freq does not matter

## k: [type = "numeric"]:
## number of genes with at least 1 count or number of positions where at least 1 read starts 

## n: [type = "numeric"]: 
## sum of counts or sum of reads mapping to the gene 
## m_boundary: [type = "numeric]:
## threshold on number size of transcript universe, e.g. number of protein coding genes

## finite sampling universe: sigma < 0 and theta = m|sigma|
## Log of Pitman is easier to compute 
 
 LogPitman_single_finite_m <- function(param,freq,k,n,m_boundary = 50000) {
  if (param[2] <= (k-1) | param[2] > m_boundary){
    return(NA)}
    else {
 
    return(sum(log(param[1] + c(1:(k-1))*(-param[1]/param[2])))-LogAscFact(param[1]+1,n -1)+sum(mapply(LogAscFact,rep(1-(-param[1]/param[2]),k),freq-1)))}
 
}


Gradient_LogPitman_single_finite_m <- function(param,freq,k,n,m_boundary = 50000) {
   if (param[2] <= (k-1)| param[2] > m_boundary){
    return(NA)}
   else {
    d <- c(0,0)
    
    d[1] <- sum(1/(param[1] + c(1:(k-1)))*(1-c(1:(k-1))/param[2])) - sum(1/(param[1]+c(1:(n-1)))) 
    	for(j in which(freq>1)) {
      		d[1] <- d[1] + sum(1/sum(param[1] + c(1:(freq[j]-1)) * param[2]))}
    
    
    d[2] <- sum(c(1:(k-1))/(param[2]^2 - c(1:(k-1))*param[2]))
    	for(j in which(freq>1)) {
      		d[2] <- d[2] + sum(-param[1]/(param[1]*param[2] + c(1:(freq[j]-1))*param[2]^2)) 
      
          }
    return(d)
  }
} 
 



# example use: optim(par=c(3000,150000),fn = LogPitman_single_infinite_m,gr  = Gradient_LogPitman_single_infinite_m, freq = freq, k = k, n = n, control = list(fnscale=-1))
 
