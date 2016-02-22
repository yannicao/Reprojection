library(geosphere)
library(RColorBrewer)

source("WRF_functions.R")

# Do some plots for the article
#

plotPolyT = function(latMin, latMax, long=0, col="red2", ...) {

  lats   = seq(latMin,latMax,.1)
  latsT  = transform2sphere(lats)
  diff   = distHaversine(cbind(long,lats),cbind(long,latsT)) / 1000 # Convert to km
  polygon(c(long,diff,long),c(latMin,lats,latMax),col=col,...)
  
}



  latMax.d03 = 42.69
  latMin.d03 = 40.97

  latMax.d02 = 44.61
  latMin.d02 = 38.93
  
  latMax.d01 = 50.23  
  latMin.d01 = 32.90
  
  lats   = 20:60
  latsT  = transform2sphere(lats)
  
  diff = distHaversine(cbind(0,lats),cbind(0,latsT)) / 1000 # Convert to km
  
  cols = brewer.pal(3,"OrRd")

  pdf("WRFLatError.pdf",width=6,height=6)
  
  plot(diff,lats,xlab="Error (km)",ylab="Latitude (deg)",pch=19,type="l")
  plotPolyT(latMin.d01, latMax.d01, col=cols[1],border=NA,density=30) 
  plotPolyT(latMin.d02, latMax.d02, col=cols[2],border=NA,density=40) 
  plotPolyT(latMin.d03, latMax.d03, col=cols[3],border=NA) 
  lines(diff,lats)
  
  legend('topright', legend=c("d01","d02","d03"),fill=cols,bty="n")
  dev.off()
  
    
# A cool website with nice great circles
#
# https://flowingdata.com/2011/05/11/how-to-map-connections-with-great-circles/