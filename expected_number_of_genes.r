# k: number of detected genes in sample of size n ( = at least one read)
# n: size of first sequencing experiment (number of reads)
# m: size of second sequencing experiment (number of reads)
# theta and sigma: parameters of the PSF estimated on the experiment of size n

expected_number_genes = function(k,n,m,theta,sigma){
i = 0:(m-1)
(prod((theta+sigma+n+i)/(theta+n+i)) -1)*(k+theta/sigma)
}

# this gives you the expected number of genes (given a first sequencing experiment with n reads and
# k detected genes) in a second sequencing experiment with s reads.