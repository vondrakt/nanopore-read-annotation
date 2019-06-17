#!/usr/bin/env Rscript

number_of_bins <- 24
size_of_bins <- 5000

satellites_to_plot <- c('rDNA_45S', 'LasTR3', 'LasTR4', 'LasTR1_6_13', 'LasTR2_16', 'LasTR5', 'LasTR7', 'LasTR8', 'LasTR9', 'LasTR10', 'LasTR11', 'LasTR12')
scatterplot_lengths <- read.table('/mnt/raid/users/tihana/190227_LAS_two_datasets/LAS180910g231_and_LAS190221/analysis_lastz_and_coding_table_190117/coded.all.scatterplot_lengths', stringsAsFactors = FALSE)

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
# collapsing lengths above the limiting bin into the last bin

for(i in 1:nrow(scatterplot_lengths)){
  
  if(scatterplot_lengths[i,2]>limiting_bin){
    scatterplot_lengths[i,2]=limiting_bin
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~
# SCATTERPLOT LENGTH PLOTTED

png('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/histogram_of_counts/Fig_histogram_of_counts_190604.png', height=2500, width=1900)

par(mar=c(6.5,5,4,6.5), oma=c(5,7.5,0,1.5), mgp=c(3,4,0))

layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), 4, 3, byrow=TRUE))

for (i in satellites_to_plot) {
  # subsetting
  scat_intact <- subset(scatterplot_lengths, V1==i & V5=='I')
  scat_truncated <- subset(scatterplot_lengths, V1==i & V5=='T')

  # histogram for intact arrays
  scat_intact_hist <- hist(scat_intact[,2], breaks = seq(from=0, to=limiting_bin, by=size_of_bins), plot=FALSE)

  # histogram for truncated arrays
  scat_truncated_hist <- hist(scat_truncated[,2], breaks = seq(from=0, to=limiting_bin, by=size_of_bins), plot=FALSE)

  barplot_data <- matrix(c(scat_intact_hist$counts, scat_truncated_hist$counts), nrow=2,byrow=TRUE)
  colnames(barplot_data) <- seq(from=0, to=limiting_bin, by=size_of_bins)[-1]
  barplot(barplot_data, col=c('red', 'blue'), names.arg = scat_intact_hist$mids, ylab='', xlab='', main=new_names[i], xaxt='n',yaxt='n',cex.main=5)

  axis(side=1, labels=c('5000', '60000', '120000'), at=c(0.6,14,28.5), line=0, lwd=5, cex.axis=5)
  axis(side=2, labels=c(pretty(0:max(barplot_data[1,]+barplot_data[,2]))), at=c(pretty(0:max(barplot_data[1,]+barplot_data[,2]))), line=0,lwd=5, cex.axis=5)

}

mtext('array length [bp]', side=1, at=0.5, line=3, outer=TRUE, cex=5)
mtext('counts', side=2, at=0.5, line=3.5,outer=TRUE,cex=5)

dev.off()
