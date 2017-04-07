  
  # Converts a shapefile to KML for use in Google Earth
  

  # Example Usage:
  # Sample_Shapefile<-readShapePoints('C:/Users/kegr/Documents/gras_projects/globwetland/toolbox/poly_smp_test_out.shp')
  # Sample_Shapefile<-readShapePoly('C:/Users/kegr/.qgis2//processing/outputs/smp_pts_out.shp')
  # crs(Sample_Shapefile)<-crs(smp)
  # 
  # Ouput_Folder<-'C:/Users/kegr/Documents/gras_projects/globwetland/toolbox/'
  # 
  # KML_Filename<-'poly_smp_test_out.kml'
  
  
  library(maptools)
  library(rgdal)
  
  vecToKml<-function(Sample_Shapefile, Ouput_Folder, KML_Filename){
  
  
    # check for .kml ectension and add if not present
    exten<-substring(KML_Filename, nchar(KML_Filename)-3, nchar(KML_Filename))
    
    if (exten !='.kml') {
      
      KML_Filename<-paste(KML_Filename, '.kml', sep="")
      
    }
    
    
    # define the output path for the kml file
    out_path<-paste(Ouput_Folder, '/' , KML_Filename, sep="")
    
    
    # reproject the shapefile to WGS84
    smp_ll<-spTransform(Sample_Shapefile, CRS("+proj=longlat +datum=WGS84"))
    
    # Get the sample lables for each point
    kml_labels<-1:nrow(Sample_Shapefile)
    
    
    # write the kml to disk
    if (class(smp_ll)[1]=='SpatialPolygonsDataFrame'){
  
    kmlPolygons(smp_ll, out_path, name=kml_label)
    #writeOGR(smp_ll,  out_path, layer='Polygon', driver="KML", overwrite_layer = TRUE)
    
    }
    
    if (class(smp_ll)[1]=='SpatialPointsDataFrame'){
      
      kmlPoints(smp_ll, out_path, name=kml_labels)
      
    }
  
  }
  
  