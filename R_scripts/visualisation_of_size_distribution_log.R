#!/usr/bin/env Rscript

# len <- read.table('/mnt/raid/users/tihana/190114_testing_LAS_data/LAS_180910/tests_db_190114/lastz_190114/pseudocoded/analysis/extracted_reads/LasTR14/LasTR14_all.sizes', stringsAsFactors = FALSE,
#                   colClasses=c("character", "numeric", "character", "numeric", "character"))
# coding_table <- read.table('/mnt/raid/454_data/vicieae/clust_jednotlive_druhy/LAS/analysis_satellites/reference_database/reference_database_190114.coding_table', stringsAsFactors = FALSE)
# limiting_bin <- 100000
# number_of_bins <- 50

args = commandArgs(trailingOnly=TRUE)

len <- read.table(args[1], stringsAsFactors = FALSE, colClasses=c("character", "numeric", "character", "numeric", "character"))
coding_table <- read.table(args[2], stringsAsFactors = FALSE)
limiting_bin <- as.numeric(args[3])
number_of_bins <- as.numeric(args[4])

pdf(args[5])

satellites <- unique(coding_table[,1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(i in 1:nrow(len)){
  
  if(len[i,2]>limiting_bin){
    len[i,2]=limiting_bin
  }
  if(len[i,4]>limiting_bin){
    len[i,4]=limiting_bin
  }
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


first_bin <- seq(from=1, to=limiting_bin, length.out = number_of_bins)[2]

# function for plotting
plot_hist <- function(x, y, start, end) {
  
  if (missing(y)) {
    
    if (x[1,5]=='I'){
      
      color='red'
      style=0
      type='intact'
      
    } else {
      
      color='blue'
      style=2
      type='truncated'
      
    }
    
    xh <- hist(x[,2], breaks = seq(from=as.numeric(start), to=as.numeric(end), length.out = number_of_bins), plot=FALSE)
    xh$counts[xh$counts==0]<-NA
    
    plot(xh$mids, xh$counts, log='x', main=paste('Size distribution of ', x[1,1], sep=''), ylab='frequency', xlab='array length', col=color, pch=style, cex=2)
    legend('topright', legend=type, col=color, pch=style)
  } else {
    
   xh <- hist(x[,2], breaks = seq(from=as.numeric(start), to=as.numeric(end), length.out = number_of_bins), plot=FALSE)
   yh <- hist(y[,2], breaks = seq(from=as.numeric(start), to=as.numeric(end), length.out = number_of_bins), plot=FALSE)
   
   xh$counts[xh$counts==0]<-NA
   yh$counts[yh$counts==0]<-NA
   
   plot(xh$mids, xh$counts, log='x', main=paste('Size distribution of ', x[1,1], sep=''), ylab='frequency', xlab='array length', col='red', pch=0, cex=2, ylim=c(1,max(c(xh$counts, yh$counts),na.rm=T)))
   points(yh$mids, yh$counts, pch=2, col='blue', cex=2)
   legend('topright', legend=c('intact', 'truncated'), col=c('red', 'blue'), pch=c(0,2), cex=2)
   
  }
  
}

# iterating over satelites
for (i in satellites) {
  
  # FIRST GRAPH
  a_first_bin <- subset(len, V1==i & V5=='I' & V2<=first_bin)
  b_first_bin <- subset(len, V1==i & V5=='T' & V2<=first_bin)
  
  if (nrow(a_first_bin)!=0 & nrow(b_first_bin)==0) {
    
    plot_hist(x=a_first_bin, start=1, end=first_bin)
    
  } else if (nrow(a_first_bin)==0 & nrow(b_first_bin)!=0) {
    
    plot_hist(x=b_first_bin, start=1, end=first_bin)
    
  } else if (nrow(a_first_bin)==0 & nrow(b_first_bin)==0) {
    
    a_rest <- subset(len, V1==i & V5=='I' & V2>first_bin)
    b_rest <- subset(len, V1==i & V5=='T' & V2>first_bin)
    
    if (nrow(a_rest)==0 & nrow(b_rest)==0) {
      
      print(paste('There are no hits for satellite ',i,sep=''))
      
    } else {
      
      print(paste('There are no hits in the first bin for satellite ',i,sep=''))
      
    }
    
  } else {
    
    plot_hist(x=a_first_bin, y=b_first_bin, start=1, end=first_bin)
    
  }
  
  # SECOND GRAPH
  
  a_rest <- subset(len, V1==i & V5=='I' & V2>first_bin)
  b_rest <- subset(len, V1==i & V5=='T' & V2>first_bin)
  
  if (nrow(a_rest)!=0 & nrow(b_rest)==0) {
    
    plot_hist(x=a_rest, start=first_bin, end=limiting_bin)
    
  } else if (nrow(a_rest)==0 & nrow(b_rest)!=0) {
    
    plot_hist(x=b_rest, start=first_bin, end=limiting_bin)
    
  } else if (nrow(a_rest)==0 & nrow(b_rest)==0) {
    
    print(paste('There are no hits outside the first bin for satellite ',i,sep=''))
    
  } else {
    
    plot_hist(x=a_rest, y=b_rest, start=first_bin, end=limiting_bin)
    
  }
  
}

dev.off()
