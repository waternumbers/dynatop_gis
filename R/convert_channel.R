#' Function for assisting in the conversion of object to be suitable channel inputs to a dynatopGIS object
#'
#' @description Converts SpatVect of the channel network into SpatVect containing polgons suitable for dynatopGIS, optionally adding in water bodies and simplifying.
#' @param chn a SpatVect object or a file which can read by terra::vect to create one
#' @param property_names a named vector of containing the columns of existing data properties required in the final SpatialPolygonsDataFrame
#' @param default_width the width in m to be used for buffering lines to produce polygons
#' @param default_slope the slope in m/m to be used when none is provided
#' @param default_depth used to compute bankfull volume is none provided
#' @param simplify_length channel reachs shorter then this will be merged.
#' @param wb a SpatVect containing waterbody polygons a or file which can read by terra::vect to create one
#'
#' @details If the property_names vector contains a width this is used for buffering lines to produce polygons, otherwise the default_width value is used.
#' wb is present must contain the feilds `name` and `length`
#' 
#' @examples
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' vect_lines <- terra::vect(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(vect_lines,property_names)
#' @export
convert_channel <- function(chn,property_names=c(name = "DRN_ID",
                                                 length = "length",
                                                 startNode = "startNode",
                                                 endNode = "endNode",
                                                 width = "width",
                                                 slope = "slope",
                                                 bankfull_vol = "volume"),
                            default_width = 2, default_slope = 0.001, default_depth = 1,
                            simplify_length = 0,
                            wb = NULL){
    
    ## read in sp object is a character sting
    if(is.character(chn)){
        if(file.exists(chn)){
            chn <- terra::vect(chn)
        }else{
            stop("chn is a character string but the file specified does not exist")
        }
    }

    if(!(terra::geomtype(chn) %in% c("lines","polygons"))){
        stop("The channel network does not have a line or polygon geometry (even when read in)")
    }
    is_polygon <- terra::geomtype(chn)=="polygons"

    ## check all feilds in property_names exist
    if( !all( property_names %in% names(chn) ) ){
        stop("A field given in property_names does not exist")
    }
    
    ## mutate the names so that they match those on the property_names
    nm <- names(chn)
    for(ii in names(property_names)){
        nm[nm == property_names[ii]] <- ii
    }
    names(chn) <- nm

    ## trim to just the required names
    chn <- chn[, names(property_names)]

    ## work out required names
    req_names <- c("name","length","startNode","endNode")
    if(!all(req_names %in% names(chn))){
        stop("Not all the required field names are present")
    }

    ## check if SpatialPolygon object - if not buffer using a width
    if(!is_polygon){
        if(!("width" %in% names(chn))){
            warning("Modifying to spatial polygons using default width")
            chn$width <- default_width
            chn <- terra::buffer(chn, width=default_width/2)
        }else{
            if(any(is.na(chn$width))){
                warning("Replacing missing channel widths with default value")
                chn$width[is.na(chn$width)] <- default_width
            }
            warning("Modifying to spatial polygons using specified width")
            chn <- terra::buffer(chn, width=chn$width/2)
        }
    }
    
    chn$area <- terra::expanse(chn)

    if(!("slope" %in% names(chn))){
        warning("Adding default slope")
        chn$slope <- default_slope
    }else{
        if(any(is.na(chn$slope))){
                warning("Replacing missing channel slopes with the default value")
                chn$slope[is.na(chn$slope)] <- default_slope
        }
    }

    if(!("bankfull_vol" %in% names(chn))){
        warning("Adding default bankfull volume")
        chn$bankfull_vol <- chn$area*default_depth
    }else{
        if(any(is.na(chn$bankfull_vol))){
                warning("Replacing missing bankfull volume with the default value")
                chn$bankfull_vol[is.na(chn$bankfull_vol)] <- chn$area[is.na(chn$bankfull_vol)] * default_depth
        }
    }

    ## some further basic checks
    chn$name <- as.character(chn$name)
    chn$startNode <- as.character(chn$startNode)
    chn$endNode <- as.character(chn$endNode)
    chn$length <- as.numeric(chn$length)
    if(!all(is.finite(c(chn$length,chn$width,chn$area,chn$bankfull_vol))) ){
        stop("Some non-finite values of channel lengths, widths, areas or volumes found!")
    }


    ## ###############################################
    ## merge in water bodies if required
    if( !is.null(wb) ){
        message("Merging waterbodies")
        if(is.character(wb)){
            if(file.exists(wb)){
                wb <- terra::vect(wb)
            }else{
                stop("wb is a character string but the file specified does not exist")
            }
        }

        stopifnot("wb must contain only spatial polygons" = terra::geomtype(wb)=="polygons",
                  "wb must have name and length fields" = all(c("name","length")%in%names(wb))
                  )
        
        wb <- wb[,c("name","length")]

        ## compute area and width and set other properties
        wb$area <- terra::expanse(chn)
        wb$width <- wb$area / wb$length
        wb$slope <- 1e-8 ## default slope
        wb$startNode <- paste0("wb_",1:nrow(wb),"_start")
        wb$endNode <- paste0("wb_",1:nrow(wb),"_end")
        wb$bankfull_vol <-  wb$area * default_depth
        wb$is_wb <- TRUE
        chn$is_wb <- FALSE

        ## intersect with channel network and replace
        for(ii in 1:nrow(wb)){
            idx <- terra::is.related(chn, wb[ii,], "intersects")
            if(!any(idx)){ next }
            ctmp <- chn[idx,]
            if(nrow(ctmp)==1){
                ## assume an outlet
                edx <- FALSE
                sdx <- TRUE
            }else{
                edx <- ctmp$endNode %in% ctmp$startNode ## true is end of reach is start of another reach in subset
                sdx <- ctmp$startNode %in% ctmp$endNode ## true is start is also and endNode for another reach in subet
            }
            
            ctmp[edx & !sdx]$endNode <- wb$startNode[ii]
            ctmp[!edx & sdx]$startNode <- wb$endNode[ii]
            ctmp <- ctmp[ !(edx & sdx) ,]
            
            chn <- chn[!idx,]
            chn <- rbind(chn,ctmp)
            chn <- rbind(chn,wb[ii,])
        }
    }

    ## ############################################################
    ## simplify the channel network
    idx <- chn$length < simplify_length
    if(any(idx)){ message("Simplifying network") }
    sN <- chn$startNode
    eN <- chn$endNode
    iidx <- which(idx)
    iidx <- iidx[order(chn$length[iidx])]
    print(length(iidx))
    cnt <- 0
    chn_df <- as.data.frame(chn)
    chn[, names(chn)] <- NULL ## reduce channel to just genometries
        
    for(ii in iidx){
        if( chn_df$length[ii] >= simplify_length ){ next }
        jj <- which( sN == eN[ii] ) ## flows out from end node
        if( length(jj) == 1 ){ ## can merge d/s since only one outflow
                                        #print(paste(ii,jj))
            ##browser()
            ##chn[[jj]] <- c( chn[[jj]], chn[[ii]])
            chn_df$length[jj] <- chn_df$length[jj] + chn_df$length[ii]
            chn_df$area[jj] <- chn_df$area[jj] + chn_df$area[ii]
            chn_df$width[jj] <- chn_df$area[jj] / chn_df$length[jj]
            chn_df$bankfull_vol[jj] <- chn_df$bankfull_vol[jj] + chn_df$bankfull_vol[ii]
            eN[eN==sN[jj]] <- sN[ii]
            sN[jj] <- sN[ii]
        }else{
            if( sum(sN == sN[ii])==1 ){ ## start node is not bifurcation so merge u/s
                jj <- which( eN==sN[ii] ) ## those terminating at start Node
                if( length(jj) > 0 ){
                    ## if length(jj)==0 there is no upstream so just remove
                    eN[jj] <- eN[ii]
                    jj <- jj[which.max(chn_df$length[jj])] ## merge into longest reach
                    #print(paste(ii,jj))
                                        #chn_df[jj,] <- terra::combineGeoms(chn_df[jj,],chn_df[ii,])
                    ##chn[[jj]] <- c( chn[[jj]], chn[[ii]])
                    chn_df$length[jj] <- chn_df$length[jj] + chn_df$length[ii]
                    chn_df$area[jj] <- chn_df$area[jj] + chn_df$area[ii]
                    chn_df$bankfull_vol[jj] <- chn_df$bankfull_vol[jj] + chn_df$bankfull_vol[ii]
                    chn_df$width[jj] <- chn_df$area[jj] / chn_df$length[jj]
                }
            }else{
                ## drop it
                #eN[ eN == sN[ii] ] <- eN[ii]
                #warning(paste("unable to merge", chn_df$name[ii], "with length", chn_df$length[ii]))
                idx[ii] <- FALSE
            }
        }
        cnt <- cnt + 1
        if( cnt/100 == cnt%/%100 ){ print(paste(cnt,"of",length(iidx))) }
    }
    ##    Reduce(combineGeoms,c(tmp[1,],tmp[2,],tmp[5,]))
    browser()
    chn <- cbind(chn,chn_df)
    chn <- chn[!idx,]

    return(chn)
}
