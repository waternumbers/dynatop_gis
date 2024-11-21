## These functions are NOT used in the core dynatopGIS functions
## They are other functions that can be generally helpfull in pre-processing

#' Locate flow gauges on a channel network
#'
#' @description Locate flow gauges on a channel network
#'
#' @param rst SpatRast (or path to file) to infill
#' @param ctch SpatRast (or path to file) defining non-NA cells
#' @param layer_name name of layer to infill
#' 
#' @details \code{rst} and \code{ctch} should have the same extent and resolution. The \code{NA} values in rst
#' are infilled using the most common neighbouring value throuhg multiple passes of focal.
#' @export
fill_na <- function(rst,ctch,layer_name = names(rst)){
    
    if(!("SpatRaster" %in% class(rst))){ layer <- terra::rast(as.character(rst)) }
    if(!("SpatRaster" %in% class(rst))){ stop("rst is not a SpatRaster") }
    if(!("SpatRaster" %in% class(ctch))){ layer <- terra::rast(as.character(ctch)) }
    if(!("SpatRaster" %in% class(ctch))){ stop("ctch is not a SpatRaster") }
    if( !all(layer_name %in% terra::names(rst)) ){ stop("layer_names are not found") }

    for(ii in layer_name){
        if( global(is.na(rst[[ii]]) & !is.na(ctch),max)>0 ){
            while( global(is.na(rst[[ii]]) & !is.na(ctch),max)>0 ){ ## if there is a big mismatch this is slow and b
                rst[[ii]] <- focal(rst[[ii]], w=3, fun="modal",na.rm=TRUE, na.policy="only",expand=TRUE)
            }
        }
    }
    rst <- mask(rst,ctch)
    return(rst)
}
 
