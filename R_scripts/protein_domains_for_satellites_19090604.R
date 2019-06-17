#!/usr/bin/env Rscript

satellites <- c('D', 'G', 'R', 'S', 'T', 'U', 'V', 'I', 'F')
domains <- c('Y','y', 'N', 'n', 'H', 'h', 'C', 'c', 'E', 'e', 'B', 'b')

codes <- c(
  'D' = 'FabTR-51',
  'G' = 'FabTR-52',
  'R' = 'FabTR-54',
  'S' = 'FabTR-55',
  'T' = 'FabTR-56',
  'U' = 'FabTR-57',
  'V' = 'FabTR-58',
  'I' = 'FabTR-59',
  'F' = 'FabTR-60',
  'Y' = 'GAG',
  'C' = 'RT',
  'H' = 'RH',
  'B' = 'aRH',
  'N' = 'INT',
  'E' = 'PROT'
)

png('/mnt/raid/users/tihana/tihana_rscipts/plots_for_manuscript/domains_for_satellite_profile/Fig_domains_for_satellite_neighborhood_profile_190604.png', height=2500, width=1900)
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9), 3,6, byrow=TRUE))
par(mar=c(6,8,7,8), oma=c(30,8,0,3), mgp=c(3,4,0))

for (i in satellites) {
  
  incrementation_table_name <- paste('/mnt/raid/users/tihana/190227_LAS_two_datasets/LAS180910g231_and_LAS190221/analysis_lastz_and_coding_table_190308/pseudocoded_Ogre_domains_190308/coded_ogre_domains.all_profiles_', i, sep='')
  base_count_table_name <- paste('/mnt/raid/users/tihana/190227_LAS_two_datasets/LAS180910g231_and_LAS190221/analysis_lastz_and_coding_table_190308/pseudocoded_Ogre_domains_190308/base_count_', i, '_coded_ogre_domains.all_profiles', sep='')
  
  # reading the tables fromfiles
  incrementation_table <- read.table(incrementation_table_name, stringsAsFactors = FALSE, header=TRUE)
  base_count_table <- read.table(base_count_table_name, stringsAsFactors = FALSE, header = TRUE)
  # removing the 0 value
  incrementation_table <- incrementation_table[,-10001]
  base_count_table <- base_count_table[,-10001]
  
  # creating the mirrored tables
  if(i != 'T') {
    
    incrementation_table <- incrementation_table[,20000:1]
    base_count_table <- base_count_table[,20000:1]
    
    color = c(
      'N' = 'red4',
      'H' = 'midnightblue',
      'B' = 'darkgoldenrod4',
      'C' = 'darkgreen',
      'E' = 'mediumpurple4',
      'Y' = 'gray',
      'n' = 'red',
      'h' = 'blue',
      'b' = 'goldenrod',
      'c' = 'green',
      'e' = 'purple',
      'y' = 'black'
      
    )
    
  } else {
    
    color = c(
      'N' = 'red',
      'H' = 'blue',
      'B' = 'goldenrod',
      'C' = 'green',
      'E' = 'purple',
      'Y' = 'black',
      'n' = 'red4',
      'h' = 'midnightblue',
      'b' = 'darkgoldenrod4',
      'c' = 'darkgreen',
      'e' = 'mediumpurple4',
      'y' = 'gray'
      )
    
  }
  
  # finding the upper limit of the y axis
  normalisation_matrix <- matrix(unlist(rep(base_count_table[1,], 12)), ncol=20000,byrow=TRUE)
  incrementation_matrix <- as.matrix(incrementation_table[c(domains),])
  ymax <- max(incrementation_matrix/normalisation_matrix, na.rm = TRUE)
  # creating the values for the x axis
  xaxis <- seq(-10000,10000)[-10001]
  
  # openning a new plot
  plot(x=xaxis, y=(incrementation_table[i,]/base_count_table[1,]), type='n', main=codes[i], ylab='', xlab='', lwd=5, bty='n', cex.main=5, xaxt='n', yaxt='n', ylim=c(0,ymax))
  axis(side=1, lwd=5, cex.axis=5, labels=c(-10000,0,10000), at=c(-10000,0,10000))
  axis(side=2, lwd=5, cex.axis=5)
  abline(v=0,lwd=5,col='black')
  
  # iterating over domains and plotting
  for (j in domains) {
    
    lines(x=xaxis, y=(incrementation_table[j,]/base_count_table[1,]), col=color[j], lwd=5)
    
  }
  
}
mtext('distance from the array [bp]', side=1, at=0.5, outer=TRUE, cex=5, line=5)
mtext('density', side=2, at=0.5, outer=TRUE, cex=5, line=1)

par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend('bottomleft', legend=c('GAG forward', 'PROT forward', 'RT forward', 'RH forward', 'aRH forward', 'INT forward',
                            'GAG reverse', 'PROT reverse', 'RT reverse', 'RH reverse', 'aRH reverse', 'INT reverse'),
       col=c('black', 'purple', 'green', 'blue','goldenrod','red','gray','mediumpurple4', 'darkgreen', 'midnightblue', 'darkgoldenrod4', 'red4'), 
       pch=c(19,19,19,19,19,19,19,19,19,19,19,19),
       cex=5,ncol=4,bty='n')


dev.off()
