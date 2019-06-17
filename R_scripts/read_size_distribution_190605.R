#!/usr/bin/env Rscript
number_of_bins <- 24
limiting_bin <- 120000
size_of_bins <- 120000/24
bins = seq(from=0, to=limiting_bin, by=size_of_bins)[-1]


png('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/read_length_distribution/Fig_read_size_distribution_190606.png', height=2500,width=1900)
layout(matrix(c(1:15), 5, 3, byrow=TRUE))
par(mar=c(7.5,6,6,6.5), oma=c(5,9,2,1.5), mgp=c(3,4,0))
options(scipen = 10)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# FIRST PLOT RUN1
run1 <- read.table('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/read_length_distribution/Run1_binned_data.table', stringsAsFactors = FALSE)

barplot(run1[,2], col='gray', names.arg = bins, ylab='', xlab='', main='Run1', cex.main=5, xaxt='n', yaxt='n')
axis(side=1, labels=c(5000,60000,120000), at=c(0.6,14,28.5), line=0, lwd=5, cex.axis=5)
numb <- pretty(0:max(run1[,2]))[-1]
axis(side=2, labels=c(0,formatC(numb,format='e',digits=2)), at=c(0,formatC(numb,format='e',digits=2)), line=0, lwd=5, cex.axis=5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SECOND PLOT RUN2
run2 <- read.table('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/read_length_distribution/Run2_binned_data.table', stringsAsFactors = FALSE)

barplot(run2[,2], col='gray', names.arg = bins, ylab='', xlab='', main='Run2', cex.main=5, xaxt='n', yaxt='n')
axis(side=1, labels=c(5000,60000,120000), at=c(0.6,14,28.5), line=0, lwd=5, cex.axis=5)
numb <- pretty(0:max(run2[,2]))[-1]
axis(side=2, labels=c(0,formatC(numb,format='e',digits=2)), at=c(0,formatC(numb,format='e',digits=2)), line=0, lwd=5, cex.axis=5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# THIRD PLOT ALL READS USED FOR ANALYSIS
all_reads <- read.table('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/read_length_distribution/All_reads_analyzed_binned_data.table', stringsAsFactors = FALSE)

barplot(all_reads[,2], col='gray', ylab='',xlab='', main='Reads used for analysis', cex.main=5, yaxt='n',xaxt='n')
axis(side=1, labels=c(5000,60000,120000), at=c(0.6,14,28.5), line=0, lwd=5, cex.axis=5)
numb <- pretty(0:max(all_reads[,2]))[-1]
axis(side=2, labels=c(0,formatC(numb,format='e',digits=2)), at=c(0,formatC(numb,format='e',digits=2)), line=0, lwd=5, cex.axis=5)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTS ITERATING OVER ALL GROUPS OF SATELLITES

satellites <- c('rDNA_45S','LasTR3', 'LasTR4', 'LasTR1_6_13', 'LasTR2_16', 'LasTR5', 'LasTR7', 'LasTR8',
                'LasTR9', 'LasTR10', 'LasTR11', 'LasTR12')

new_names <- c(
  'rDNA_45S' = '45S rDNA',
  'LasTR4' = 'FabTR-53',
  'LasTR3' = 'FabTR-2',
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

for(i in satellites) {
  
  file_name = paste('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/read_length_distribution/Satellite_binned_data.table_',i,sep='')
  read_group <- read.table(file_name, stringsAsFactors = FALSE)
  
  barplot(read_group[,2], col='gray', names.arg = bins, ylab='', xlab='', main=new_names[i], cex.main=5, xaxt='n', yaxt='n')    
  axis(side=1, labels=c(5000,60000,120000), at=c(0.6,14,28.5), line=0, lwd=5, cex.axis=5)
  numb <- pretty(0:max(read_group[,2]))[-1]
  axis(side=2, labels=c(0,formatC(numb,format='e',digits=2)), at=c(0,formatC(numb,format='e',digits=2)), line=0, lwd=5, cex.axis=5)
  
}

mtext('read length [bp]', side=1, at=0.5, line=3, outer=TRUE, cex=5)
mtext('sum of length [bp]', side=2, at=0.5, line=4,outer=TRUE,cex=5)

dev.off()
