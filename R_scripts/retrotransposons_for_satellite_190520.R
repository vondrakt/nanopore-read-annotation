#!/usr/bin/env Rscript

satellites <- c('A', 'O', 'P', 'D', 'G', 'R', 'S', 'T', 'U', 'V', 'I', 'F')
retrotransposons <- c('C', 'E', 'Y', 'Z')

color = c(
  'B' = 'sienna2',
  'C' = 'green',
  'W' = 'goldenrod',
  'Q' = 'purple',
  'E' = 'plum1',
  'Y' = 'blue',
  'H' = 'salmon4',
  'Z' = 'darkgreen',
  'r' = 'red',
  't' = 'plum4'
)

codes <- c(
  'A' = '45S rDNA',
  'P' = 'FabTR-53',
  'O' = 'FabTR-2',
  'D' = 'FabTR-51',
  'G' = 'FabTR-52',
  'R' = 'FabTR-54',
  'S' = 'FabTR-55',
  'T' = 'FabTR-56',
  'U' = 'FabTR-57',
  'V' = 'FabTR-58',
  'I' = 'FabTR-59',
  'F' = 'FabTR-60',
  'B' = 'DNA',
  'C' = 'LTR_Copia_Maximus',
  'W' = 'LTR_Copia_other',
  'Q' = 'LTR_gypsy_Athila',
  'E' = 'LTR_gypsy_chromo',
  'Y' = 'LTR_gypsy_Ogre',
  'H' = 'LTR_gypsy_other',
  'Z' = 'LTR_other'
)


png('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/retrotransposon_to_satellite_profile/Fig_retrotransposon_to_satellite_profile_190612.png',  height=2500, width=1900)

layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12), 4,6, byrow=TRUE))

par(mar=c(6,8,7,8), oma=c(23,8,0,1.5), mgp=c(3,4,0))

for (i in satellites) {
  
  incrementation_table_name <- paste('/mnt/raid/users/tihana/190227_LAS_two_datasets/LAS180910g231_and_LAS190221/analysis_lastz_and_coding_table_190306/coded.all.profiles_', i, sep='')
  base_count_table_name <- paste('/mnt/raid/users/tihana/190227_LAS_two_datasets/LAS180910g231_and_LAS190221/analysis_lastz_and_coding_table_190306/base_count_', i, '_coded.all.profiles', sep='')
  # loading the table from file
  incrementation_table <- read.table(incrementation_table_name, stringsAsFactors = FALSE, header=TRUE)
  base_count_table <- read.table(base_count_table_name, stringsAsFactors = FALSE, header = TRUE)
  # removing the 0 position from the table that caused a break in the line graph
  incrementation_table <- incrementation_table[,-10001]
  base_count_table <- base_count_table[,-10001]
  # creating a mirrored table 
  if (i != 'T') {
  
    incrementation_table<-incrementation_table[,20000:1]
    base_count_table <- base_count_table[,20000:1]
  
  }
  # creating the x axis without 0
  xaxis <- seq(-10000,10000)[-10001]
  
  # plotting the satellite to itself

  plot(x=xaxis, y=(incrementation_table[i,]/base_count_table[1,]), type='l', col= 'black', main=codes[i], ylab='', xlab='', lwd=1, bty='n', cex.main=1, xaxt='n', yaxt='n', ylim=c(0,1), cex.main=5)
  lines(x=xaxis, y=(incrementation_table[tolower(i),]/base_count_table[1,]), col= 'gray', lwd=5)
  axis(side=1, lwd=5, cex.axis=5, labels=c(-10000,0,10000), at=c(-10000,0,10000))
  axis(side=2, lwd=5, cex.axis=5)
  if (any(incrementation_table[i,]/base_count_table[1,] > 0.1, na.rm=TRUE)) lines(x=xaxis, y=(incrementation_table[i,]/base_count_table[1,]), col='black', lwd=5)
  if (any(incrementation_table[tolower(i),]/base_count_table[1,] > 0.1, na.rm=TRUE)) lines(x=xaxis, y=(incrementation_table[tolower(i),]/base_count_table[1,]), col= 'gray', lwd=5)
  abline(v=0,lwd=5,col='black')
  lines(x=xaxis, y=(incrementation_table['r',]/base_count_table[1,]), col= color['r'], lwd=5)
  
  for (j in retrotransposons) {
    
    lines(x=xaxis, y=(incrementation_table[j,]/base_count_table[1,]), col=color[j], lwd=5)
    
  }
  
}
mtext('distance from the array [bp]', side=1, at=0.5, outer=TRUE, cex=5, line=5)
mtext('density', side=2, at=0.5, outer=TRUE, cex=5, line=2)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')

legend('bottom', legend=c('same satellite as examined', 'Ty3/gypsy/Ogre', 'LTR/unclassified', 'Ty3/gypsy/Chromovirus','Ty1/copia/Maximus', 'FabTR-54'),
       col=c('black','blue', 'darkgreen', 'plum2', 'green','red'),
       pch=c(15,15,15,15,15),
       cex=5, ncol=3,bty='n')

dev.off()
