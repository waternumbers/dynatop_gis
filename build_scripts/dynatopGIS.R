rm(list=ls())
## my_new_object  <-  my_object

## my_new_object  <-  my_object$clone()

devtools::load_all()
#library("dynatopGIS")

demo_dir <- "./build_scripts/demo" #tempfile("dygis")
unlink(demo_dir,recursive=TRUE)
dir.create(demo_dir)

ctch <- dynatopGIS$new(file.path(demo_dir))

dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopGIS", mustWork = TRUE)

dem <- terra::rast(dem_file)
dem <- terra::extend(dem,1) ## pad with NA values
catchment_outline <- terra::ifel(is.finite(dem),1,NA)
ctch$add_catchment(catchment_outline)

ctch$add_dem(dem)

sp_lines <- terra::vect(channel_file)
head(sp_lines)

property_names <- c(name="identifier",
                    endNode="endNode",
                    startNode="startNode",
                    length="length")
chn <- convert_channel(sp_lines,property_names)

ctch$add_channel(chn)

ctch$get_layer()

ctch$plot_layer("dem", add_channel=TRUE)

ctch$get_layer("dem")

ctch$sink_fill()

terra::plot( ctch$get_layer('filled_dem') - ctch$get_layer('dem'),
            main="Changes to height")

ctch$compute_band()
ctch$plot_layer("band")

tmp <- ctch$get_layer("filled_dem")
tmp <- terra::ifel(tmp<=500,NA,-999)

ctch$add_layer(tmp, "greater_500")
ctch$get_layer()


ctch$create_model(file.path(demo_dir,"new_model"),verbose=TRUE)

list.files(demo_dir,pattern="new_model*")

#tmp <- readRDS( file.path(demo_dir, "new_model.rds") )
#mdl <- dynatop::dynatop$new(tmp$hru)

## ## ## ###################################
## rm(list=ls())
## mdl <- readRDS("./demo/new_model.rds")$hru
## bnd <- sapply(mdl,function(x){x$band})
## for(h in mdl){
##     if( length(h$sf_flow_direction$id) != length(h$sf_flow_direction$fraction) ){
##         browser()
##     }
## }



## ## ###################################
## rm(list=ls())
## mdl <- readRDS("./demo/new_model.rds")$hru
## bnd <- sapply(mdl,function(x){x$band})
## for(h in mdl){
##     if(any( bnd[h$sf_flow_direction$id+1] >= h$band )){
##         browser()
##     }
## }


## ## ###################################
## rm(list=ls())
## library(terra)

## bnd <- rast( "./demo/band.tif")

## mbnd <- as.matrix(bnd,wide=TRUE)

## dr <- readRDS("./demo/dem.rds")

## nr <- nrow(mbnd)
## delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

## for(nr in nrow(dr)){
##     ii <- dr[nr,1] ##cell
##     jj <- ii+delta
##     jj <- jj[ dr[nr,4:11]>0 ]

##     if( !all( mbnd[jj] < mbnd[ii] ) ){ stop("Oh no!!") }
## }
