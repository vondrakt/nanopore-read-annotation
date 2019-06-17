#!/usr/bin/env Rscript

# This script takes the tabular binning data from python script:
# /mnt/raid/users/tihana/tihana_pyscripts/180813/plotting_cumulative_lengths_and_frequency_of_occurences/plotting_cumulative_lengths_and_frequency_of_occurences_02.py
# and plots the weighted histogram.

# parsing arguments passed to the script
args = commandArgs(trailingOnly=TRUE)

path=args[1]
pattern=args[2]

# files_of_length will provide the patch to the output and pattern
# according to witch the output of the python script will be found.
files_length = list.files(path = path, pattern = pattern)
size_of_bins = as.numeric(args[3])
number_of_bins = as.numeric(args[4])

pdf(args[5])

# calculating the limiting bin
limiting_bin <- size_of_bins * number_of_bins
# getting the bin ranges
bins <- seq(from=size_of_bins, to=limiting_bin, by=size_of_bins)

# iterating over files of binning data for all groups of satellites
for(i in files_length) {

  # reading the tabular input separately for each satellite
  length_table <- read.table(paste(args[1],'/',i,sep=''),stringsAsFactors = FALSE, header=TRUE)
  # creating the barplot
  barplot(t(length_table[,1:2]), ylim= c(0,max(length_table[,1]+max(length_table[,2]))), col=c('red', 'blue'), main=i, names.arg = bins, las=2)
  # creating a legend
  legend('top', legend=c('inact', 'truncated'), pch=c(15,15), col=c('red', 'blue'))
  
  
}

dev.off()
