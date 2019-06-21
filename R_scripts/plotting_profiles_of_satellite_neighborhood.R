#!/usr/bin/env Rscript

# This script takes the tabular output of the python script:
# /mnt/raid/users/tihana/tihana_pyscripts/180813/profiles_of_neighborhoods/profile_of_neighborhood_04.py
# and plots the profiles of neighborhoods.

args = commandArgs(trailingOnly=TRUE)

pat = paste(args[2],'_.*', sep = '')
# files_distribution and files_counted_bases will provide the path to the output and pattern
# according to witch the output of the python script will be found.
files_distribution = list.files(path = args[1], pattern = pat)
files_counted_bases = list.files(path = args[1], pattern = args[3])

pdf(args[4])

window_size = as.integer(args[5])

# plotting the profiles of the neighborhood

# dividing the plotting area into four individual plots
par(mfrow = c(2,2))

# iterating over files
for (i in files_distribution){

  file_name <- paste(args[1], '/', i, sep = '')
  # reading the tabular input
  table <- read.table(file_name, stringsAsFactors = FALSE, header = TRUE)
  # row_names stores the names of the pseudocodes for whichthe neighborhood is to be plotted
  row_names <- rownames(table)

  # iterating over all rows of a table
  for (j in 1:nrow(table)){
    
    # if the row of the table is not filled with only 0
    if(any(table[j,]>0)) {
      # create the plot name
      plot_name = paste('Profile of ',row_names[j], ' for ',i)
      # create the linear plot of the neighborhood
      plot(x=seq(-window_size,window_size), y=table[j,], type='l', xlab='positions', ylab=row_names[j], main=plot_name, cex.main = 0.5)
    }
    
  }

}

# plotting the distribution of plotted bases

# iterating over files
for (i in files_counted_bases) {

  file_name = paste(args[1], '/', i, sep = '')
  # reading the tabular input
  table = read.table(file_name, stringsAsFactors = FALSE, header = TRUE)
  # create the plot name
  plot_name = paste('Profile of plotted bases for ',i, sep = '')

  # iterating over all rows of a table
  for (j in 1:nrow(table)) {

    # create the linear plot of distribution
    plot(x=seq(-window_size,window_size), y=table[j,], type = 'l', xlab='positions', ylab='number of plotted bases', main=plot_name, cex.main = 0.5)

  }

}

# writting to pdf

dev.off()
