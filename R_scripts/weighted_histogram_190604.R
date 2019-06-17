#!/usr/bin/env Rscript

number_of_bins <- 24
size_of_bins <- 5000

satellites_to_plot <- c('rDNA_45S', 'LasTR3', 'LasTR4', 'LasTR1_6_13', 'LasTR2_16', 'LasTR5', 'LasTR7', 'LasTR8', 'LasTR9', 'LasTR10', 'LasTR11', 'LasTR12')

new_names <- c(
  'rDNA_45S' = '45S rDNA',
  'LasTR3' = 'FabTR-2',
  'LasTR4' = 'FabTR-53',
  'LasTR1_6_13' = 'FabTR-51',
  'LasTR2_16' = 'FabTR-52',
  'LasTR5' = 'FabTR-54',
  'LasTR7' = 'FabTR-55',
  'LasTR8' = 'FabTR-56',
  'LasTR9' = 'FabTR-57',
  'LasTR10' = 'FabTR-58',
  'LasTR11' = 'FabTR-59',
  'LasTR12' = 'FabTR-60'
)

limiting_bin <- size_of_bins*number_of_bins
#~~~~~~~~~~~~~~~~~~~~~~
# LAYOUT OF THE PLOT

png('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/weighted_length_histogram/Fig_weighted_length_histograms_190604.png', height=2500, width=1900)

# par(mar=c(4,3,2,0), cex.axis=5, cex.lab=5, cex.main=5,oma=c(3,3,0,1))
par(mar=c(7.5,5,6,6.5), oma=c(5,9,2,1.5), mgp=c(3,4,0))

layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), 4, 3, byrow=TRUE))

# PLOTTING THE LENGTH WEIGHTED HISTOGRAM

for (i in satellites_to_plot) {
  file_name <- paste('/mnt/raid/users/tihana/190227_LAS_two_datasets/LAS180910g231_and_LAS190221/analysis_lastz_and_coding_table_190117/coded.all.cumulative_length_', i, sep='')
  cumulative_lengths <- read.table(file_name, stringsAsFactors = FALSE, header=TRUE)
  
  bins = seq(from=0, to=limiting_bin, by=size_of_bins)[-1]
  options(scipen=0)
  barplot_data <- t(cumulative_lengths[,1:2])
  colnames(barplot_data) <- seq(from=0, to=limiting_bin, by=size_of_bins)[-1]
  barplot(barplot_data, col=c('red', 'blue'), names.arg = bins, ylab='', xlab='', main=new_names[i], xaxt='n', yaxt='n', cex.main=5)
  
  axis(side=1, labels=c('5000', '60000', '120000'), at=c(0.6,14,28.5), line=0, lwd=5, cex.axis=5)
  numb <- pretty(0:max(cumulative_lengths[,1]+cumulative_lengths[,2]))[-1]
  axis(side=2, labels=c(0,formatC(numb,format='e',digits=2)), at=c(0,formatC(numb,format='e',digits=2)), line=0, lwd=5, cex.axis=5)
}

mtext('array length [bp]', side=1, at=0.5, line=3, outer=TRUE, cex=5)
mtext('sum of length [bp]', side=2, at=0.5, line=4,outer=TRUE,cex=5)

dev.off()
