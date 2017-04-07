##[Sampling]=group
##Sample_Shapefile=vector
##Ouput_Folder=Folder
##KML_Filename=string enter_file_name_here.kml

#tryCatch(find.package("maptools"), error=function(e) install.packages("maptools", lib=file.path(.Library[1])))

library(maptools)

Rpath <- file.path(Sys.getenv("USERPROFILE"), ".qgis2", "processing", "rscripts", fsep="\\")

fun_path <-paste(Rpath, '\\vecToKml.R', sep='')

source(fun_path)

vecToKml(Sample_Shapefile, Ouput_Folder, KML_Filename)
