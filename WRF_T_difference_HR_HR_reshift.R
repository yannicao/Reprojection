
#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
# File: WRF_functions.R
#
#Author: Guido Cervone (cervone@psu.edu) and Yanni Cao (yvc5268@psu.edu)
#        Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#        Department of Geography and Institute for CyberScience
#        The Pennsylvania State University
#
#



# This file is used to write the output difference once WRF has completed the runs

require(rgdal)
library(ncdf)
require(tools)
require(raster)

#using the new source code "WRF_function_result_Yanni"since WRF output files uses XLAT, XLAT
#describe CLAT,CLONG in WRF input file
#
source("~/geolab/projects/Marcellus/Reprojection/WRF_functions.R")
source("~/geolab/projects/Rbin/RGIS_functions.R")

CRS.WGS84  = CRS("+proj=longlat +datum=WGS84")


getDataWRF = function( filename.wrf, WRF.layer="Temp", layer=2 ) {
  
  print(paste("Processing:",filename.wrf,WRF.layer,layer))
  
  WRF                    = open.ncdf(filename.wrf,write=F)

  T               = get.var.ncdf(WRF,"T")[,,layer]
  
  if ( WRF.layer == "Temp") {
    #Calculate temprature, combing:
    #http://gradsusr.org/pipermail/gradsusr/2011-December/031698.html 
    ###this one is wrong since WRF Pressure is in pa
    ###it should be converted to mbar. so I used (PB+B)/100
    #http://mailman.ucar.edu/pipermail/wrf-users/2010/001896.html
    P               = get.var.ncdf(WRF,"P")[,,layer]
    PB              = get.var.ncdf(WRF,"PB")[,,layer]
    TotalPotentialTemperature = T+300
    WRF.data     = TotalPotentialTemperature*((((PB+P)/100)/1000)^(2/7))-273.15
  } else {
    WRF.data     = get.var.ncdf(WRF,WRF.layer)[,,layer]
    WRF.data     = WRF.data[1:nrow(T),1:ncol(T)]
  }
  

  filename.wrf.grid      = paste(filename.wrf,"_grid.Rdata",sep="")
  
  if ( !file.exists(filename.wrf.grid) ) {
    WRF.grid.polygons  = WRF.grid2polygons(WRF, type="output")
    WRF.grid.sp        = SpatialPolygons(WRF.grid.polygons, proj4string = CRS.WGS84)
    save(WRF.grid.sp, file=filename.wrf.grid)
  } else {
    load(filename.wrf.grid)
  }
  
  res.df                 = data.frame(CellID = 1:length(WRF.data), Data=as.vector(t(WRF.data)))   
  WRF.grid.spdf          = SpatialPolygonsDataFrame(WRF.grid.sp, data=res.df)
  
  return( list(WRF.grid.spdf, WRF.data ) )
}


getRasterWRF = function( filename.wrf, WRF.layer="Temp", layer=2 ) {

  temp           = getDataWRF(filename.wrf, WRF.layer, layer)
  WRF.grid.spdf  = temp[[1]]
  WRF.data       = temp[[2]]
  
  raster.temp            = flip(raster(t(WRF.data)),direction="y")
  extent(raster.temp)    = extent(WRF.grid.spdf)
  ncol(raster.temp)      = nrow(WRF.data)  # thins in R are rotated and flipped compared to raster
  nrow(raster.temp)      = ncol(WRF.data)
  
  #return(raster.temp)
  
  data.raster            = rasterize( WRF.grid.spdf, raster.temp, field='Data', fun=mean)
  
  #writeRaster(data.raster, filename = "temp.tif", format="GTiff")
  #CRS.WRF    = CRS("+proj=lcc +lat_1=30 +lat_2=60 +lat_0=41.8389129639 +lon_0=-77 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314 +units=m +no_defs")
  #gdalwarp("temp.tif", "tempout.tif",t_srs=(CRS.WGS84), t_srs=(CRS.WRF), verbose=TRUE, overwrite = FALSE) 
  #tempout    = raster("tempout.tif")
  #extent(tempout) = extent(WRF.grid.spdf)
  
  return (data.raster)
}



processWRF = function(dir1, dir2, dir3, id, WRF.layer, layer ) {

  print(paste(dir1, dir2, dir3, id, WRF.layer, layer))
  
  files                  = list.files(path=dir1, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR              = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  files                  = list.files(path=dir2, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR.SHIFT        = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  files                  = list.files(path=dir3, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.DEFAULT         = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  # Shift up (reshift)
  #
  raster.HR.RESHIFT      = raster.HR.SHIFT
  extent(raster.HR.RESHIFT)[3:4] = transform2spheroid( extent(raster.HR.RESHIFT)[3:4])
  
  
  extent1 = extent(raster.HR)
  extent2 = extent(raster.HR.RESHIFT)
  
  extent3 = extent(extent1[1], 
                   extent1[2], 
                   max(extent1[3],extent2[3]), 
                   min(extent1[4],extent2[4]) )
  spoly   = extent2SpatialPolygons( extent3, CRS.WGS84 )
  
  raster.HR.RESHIFT.crop = crop(raster.HR.RESHIFT, spoly, snap="out" )
  raster.HR.SHIFT.crop   = crop(raster.HR.SHIFT, spoly, snap="out" )
  raster.HR.crop         = crop(raster.HR, spoly, snap="out" )
  raster.DEFAULT.crop    = crop(raster.DEFAULT, spoly, snap="out" )
  
  diff1 = raster.HR.crop - raster.DEFAULT.crop
  diff2 = raster.HR.crop - raster.HR.SHIFT.crop
  diff3 = raster.HR.crop - raster.HR.RESHIFT.crop
  
  
  writeRaster(raster.HR.crop, file=paste("raster.HR.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(raster.HR.SHIFT.crop, file=paste("raster.HR.SHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(raster.HR.RESHIFT.crop, file=paste("raster.HR.RESHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(raster.DEFAULT.crop, file=paste("raster.DEFAULT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  
  writeRaster(diff1, file=paste("diff.HR-DEFAULT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(diff2, file=paste("diff.HR-HR.SHIFT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(diff3, file=paste("diff.HR-HR.RESHIFT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
}



processWRFFast = function(dir1, dir2, id, WRF.layer, layer ) {
  
  print(paste(dir1, dir2, dir3, id, WRF.layer, layer))
  
  files                  = list.files(path=dir1, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR              = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  files                  = list.files(path=dir2, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR.SHIFT        = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  # Shift up (reshift)
  #
  raster.HR.RESHIFT      = raster.HR.SHIFT
  extent(raster.HR.RESHIFT)[3:4] = transform2spheroid( extent(raster.HR.RESHIFT)[3:4])
  
  
  extent1 = extent(raster.HR)
  extent2 = extent(raster.HR.RESHIFT)
  
  extent3 = extent(extent1[1], 
                   extent1[2], 
                   max(extent1[3],extent2[3]), 
                   min(extent1[4],extent2[4]) )
  spoly   = extent2SpatialPolygons( extent3, CRS.WGS84 )
  
  raster.HR.crop         = crop(raster.HR, spoly, snap="out" )
  raster.HR.RESHIFT.crop = crop(raster.HR.RESHIFT, spoly, snap="out" )
  diff3 = raster.HR.crop - raster.HR.RESHIFT.crop
  
  writeRaster(raster.HR.RESHIFT.crop, file=paste("raster.HR.RESHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)

  writeRaster(diff3, file=paste("diff.HR-HR.RESHIFT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
}


#dir1="~/geolab/data/Marcellus/Reprojection/WRF_result/HR"
dir1="~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Feb2016/d03_HR/"
#dir2="~/geolab/data/Marcellus/Reprojection/WRF_result/HR_shift"
dir2="~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Feb2016/d03_HR_shift/"
dir3="~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Dec2015/Default/"


layers     = c(2,7,13, 19, 25)     # The hour for the result (e.g. 30 = 30th hour from beginning of simulation)
WRF.layers = c("Temp","U","V")
ids        = c(1,13,25)          # vertical layer for the parameter (2 = second layer from surface)

for ( layer in layers ) {
  for (WRF.layer in WRF.layers) {
    for (id in ids) {
      processWRFFast(dir1, dir2,  id, WRF.layer, layer)      
    }
  }
}




