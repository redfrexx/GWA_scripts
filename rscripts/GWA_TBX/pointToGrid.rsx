##[Sampling]=group
##Point_Samples=vector
##Raster_To_Base_Grid_On=raster

##Sample_Outline=output vector

Rpath <- file.path(Sys.getenv("USERPROFILE"), ".qgis2", "processing", "rscripts", fsep="\\")

fun_path <-paste(Rpath, '\\pointToGrid.R', sep='')

source(fun_path)

Result <- PointToGrid(Raster_To_Base_Grid_On, Point_Samples)

Sample_Outline<-Result
