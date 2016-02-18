library(raster)
library(stringi)

dir="/Volumes/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Comparisons-Feb2016"
setwd(dir)

str_name <- list.files(pattern = "diff.HR-HR.RESHIFT.Temp_.*\\.tif", recursive=TRUE)

raster <- raster(str_name[3])
n<- na.omit(getValues(raster))
col.name = stri_sub(str_name[3],25,-5)
dat = data.frame(data=n, names = c(rep(col.name, length(n))))

order = c(6,11,2,7,12,3,8,13,4,9,14,5,10,15)
order = c(8,13,5,10,15,1,6,11,2,7,12,4,9,14)
for (i in order) {
  str_name[i]
  raster <- raster(str_name[i])
  n<- na.omit(getValues(raster))
  col.name = stri_sub(str_name[i],25,-5)
  dat = rbind(dat,data.frame(data=n,names=col.name))
  #hist(n,main = i,breaks="FD")
}

boxplot(data ~ names, data =dat, 
        outline = FALSE,
        ylab ="Temprature Difference HR VS HR_Reshift", 
        xlab ="",
        las =2, 
        col = c("red","sienna","palevioletred1"),
        at = c(1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19),
        names = c("May 14th 12PM", "May 14th 23PM", "May 15th 12PM","May 14th 12PM", "May 14th 23PM", "May 15th 12PM","May 14th 12PM", "May 14th 23PM", "May 15th 12PM","May 14th 12PM", "May 14th 23PM", "May 15th 12PM", "May 14th 12PM", "May 14th 23PM", "May 15th 12PM"),
        par(mar = c(12, 5, 4, 2)+ 0.1)
        )
mtext("Difference",1,line=8)

