# k: number of detected genes in present sample of size n
# n: size of first sequencing sample (number of reads)
# m: size of second sequencing sample (number of reads)
# theta and sigma: parameters of the PSF calibrated on the sample of size n

expected_number_genes = function(k,n,m,theta,sigma){
i = 0:(m-1)
(prod((theta+sigma+n+i)/(theta+n+i)) -1)*(k+theta/sigma)
}
