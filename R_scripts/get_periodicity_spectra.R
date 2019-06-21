#!/usr/bin/env Rscript
library(Biostrings)
library(TSclust)
library(Rfast)
convert.fft <- function(cs, sample.rate=1) {
  cs <- cs / length(cs) # normalize

  distance.center <- function(c)signif( Mod(c),        4)
  angle           <- function(c)signif( 180*Arg(c)/pi, 3)
  df <- data.frame(cycle    = 0:(length(cs)-1),
                   freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                   strength = sapply(cs, distance.center),
                   delay    = sapply(cs, angle))
  df$monomer = nrow(df)/df$cycle
  df
}

dna2vector=function(s, method="basic"){
  ## assume single sequence!
  coding1 = c(A = -1, C = 1, G = -2 ,T = -2)
  coding2 = c(A = -1, C = 0, G = 0 ,T = 1)
  coding3 = c(A = 0, C = -1, G = 1 ,T = 0)

  if (method == "basic"){
    x = as.numeric(as.factor(unlist(strsplit(as.character(s),split=""))))
  }
  if (method == "pupy"){
    x = coding1[unlist(strsplit(as.character(s),split=""))]
  }
  if (method == "cumulative"){
    x = cumsum(coding1_corrected[unlist(strsplit(as.character(s),split=""))])
    ## detrend:
    x = x- lowess(x, f= 0.01)$y
  }
  if (method=="pupy-complement"){
    x1 = coding2[unlist(strsplit(as.character(s),split=""))]
    x2 = coding3[unlist(strsplit(as.character(s),split=""))]
    x = cbind(x1,x2)
  }
  return(x)
}

get_periodicity = function(x, window=5000, step=10, xrange=seq(1,15000,1)){
  L=length(x)
  dflist=list()
  index=numeric()
  j=0
  y = numeric(length(xrange))
  ## for (i in seq(0,window,step)){
  for (i in unique((round(1.4^(seq(0,log(window,base=1.4))))))){
    print(i)
    j = j + 1
    df = convert.fft(fft(x[1:(L-i)]))
    y = y + approx(x = df$monomer[-1], y=df$strength[-1], xout = xrange)$y
  }
  y/j
}


ACF2 = function(s,L){
  x = unlist(strsplit(tolower(unname(as.character(s))),""))
  S1_4=list()
  for (j in c("a","c","g","t")){
    A = as.numeric(x == j)
    S1_4[[j]] = acf(A,lag.max = L, plot=FALSE)$acf[-1,1,1]
  }
  do.call(cbind, S1_4)
}

## input fasta file
sfile = commandArgs(T)
print(sfile)
q()

## get periodicity using fft
for (f in  sfile){
  print(f)
  s = readDNAStringSet(f)
  sbig = s[nchar(s)>=30000]
  xbig = lapply(sbig, dna2vector, method = "basic")
  periodicity = mclapply(xbig, FUN = get_periodicity, mc.cores = 40, window=200)
  periodicity_matrix = do.call(cbind,periodicity )
  saveRDS(periodicity_matrix, file=paste0(f,"_periodicity_fft.RDS"))
}

## get periodicity using autocorrelation
for (f in  sfile){
  print(f)
  s = readDNAStringSet(f)
  sbig = s[nchar(s)>=30000]
  periodicity = lapply(sbig, FUN = ACF2, L=10000)
  saveRDS(periodicity, file=paste0(f,"_periodicity_acf.RDS"))
}

