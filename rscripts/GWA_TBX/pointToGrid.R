  
  #---------------------------------------------------------------------------------#

  # This function takes a point shapefile and makes a square grid for every point
  # based on the grid of a specified raster.

  # For each point that fall within a raster pixel cell, a square ploygon is made 
  # with the outline being the raster pixel

  #---------------------------------------------------------------------------------#

  
  # function to convert a pixel (cell) to polygon
  # can be used in a loop if you have a vector of cell numbers
  # x = base raster
  # cell = cell number of the base raster (these are converted to polygons)  
  PixToPoly<-function (cell, x) {
    
    # get the xy position
    cell <- as.vector(xyFromCell(x, cell))
    
    # get the center offset
    off<-(res(x)/2)[1]
    
    # make the polygon - accounting for the center offset
    e <- extent(c(cell[1]-off, cell[1]-off + res(x)[1], cell[2]-off, cell[2]-off + res(x)[2]))
    e <- as(e, "SpatialPolygons")
    projection(e) <- projection(x)
    return(e)
  }
  
  
  PointToGrid<-function(inRst, smp){
    
    # get attribute data
    att<-smp@data
    
    # get the cell number raster corresponding to each points
    clnum<-as.matrix(cellFromXY(inRst, smp))
    
    # apply PixToPoly function to create grid polygons based on cell numbers
    pol_lst<-apply(clnum,  1, PixToPoly, x=inRst)
    
    # merge the list of grid polygons
    pol<-do.call(bind, pol_lst)
    
    # add attribute data to pixel polygons  
    pol_spd<-SpatialPolygonsDataFrame(pol, data=att)
    
    return(pol_spd)
    
  }
  