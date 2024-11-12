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
convert_channel <- function(chn,property_names=c(name = "name",
                                                 length = "length",
                                                 startNode = "startNode",
                                                 endNode = "endNode",
                                                 width = "width",
                                                 slope = "slope",
                                                 channelVol = "channelVol"),
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

    if(!("channelVol" %in% names(chn))){
        warning("Adding default bankfull volume")
        chn$channelVol <- chn$area*default_depth
    }else{
        if(any(is.na(chn$channelVol))){
                warning("Replacing missing bankfull volume with the default value")
                chn$channelVol[is.na(chn$channelVol)] <- chn$area[is.na(chn$channelVol)] * default_depth
        }
    }

    ## some further basic checks
    chn$name <- as.character(chn$name)
    chn$startNode <- as.character(chn$startNode)
    chn$endNode <- as.character(chn$endNode)
    chn$length <- as.numeric(chn$length)
    if(!all(is.finite(c(chn$length,chn$width,chn$area,chn$channelVol))) ){
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
        wb$channelVol <-  wb$area * default_depth
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
    if( any(chn$length < simplify_length) ){
        message("Simplifying network")
        gm <- chn
        gm[] <- NULL
        gm <- split(gm,1:nrow(gm))
        chn <- as.data.frame(chn)
        to_keep <- rep(TRUE,nrow(chn))
        
        ## find channels to merge and order
        idx <- which( chn$length < simplify_length )
        idx <- idx[order(chn$length[idx])]

        cnt <- 0
        n <- length(idx)
        print(paste("Evaluating",n,"potential candidates"))
        pb <- txtProgressBar(min = 0, max = n, initial = 0, char = "=",
                             width = NA, "progress...", "whoop", style = 3, file = "")
        for(ii in idx){

            cnt <- cnt+1
            setTxtProgressBar(pb, cnt, title = NULL, label = NULL)
            
            if( chn$length[ii] >= simplify_length ){ next } ## catch incase a merge has already made it long
            
            in_hn <-  which( chn$endNode == chn$startNode[ii] )
            out_hn <-  which( chn$startNode == chn$startNode[ii] )
            in_en <- which( chn$endNode == chn$endNode[ii] )
            out_en <- which( chn$startNode == chn$endNode[ii] )
            jj <- NA
            if( length(in_en)==1 &
                length(out_en)==1 ){
                ## endNode is a 1-to-1 join and can merge d/s
                jj <- out_en
                chn$startNode[jj] <- chn$startNode[ii]
            }
            if( length(in_hn)==1 &
                length(out_hn)==1 &
                is.na(jj) ){
                ## startNode is a 1-to-1 join and can merge u/s
                jj <- in_hn
                chn$endNode[jj] <- chn$endNode[ii]
            }
            if( length(out_en)==1 &
                is.na(jj) ){
                ## endNode is a single outlet so merge d/s but reroute other flows to endNode
                jj <- out_en
                chn$endNode[ chn$endNode==chn$startNode[jj] ] <- chn$startNode[ii]
                chn$startNode[jj] <- chn$startNode[ii]
            }
            
            if(is.na(jj)){ next } ## can't simplify

            ##print(paste(cnt,ii,jj))
            
            ## merge properties (default to those for jj)
            chn$length[jj] <- chn$length[jj] + chn$length[ii]
            chn$area[jj] <- chn$area[jj] + chn$area[ii]
            chn$width[jj] <- chn$area[jj] / chn$length[jj]
            chn$channelVol[jj] <- chn$channelVol[jj] + chn$channelVol[ii]
            ## chn[ii,] <- NA # usefull for debugging...

            ## merge geom
            gm[[jj]] <- combineGeoms( gm[[jj]], gm[[ii]], minover=0 )
            to_keep[ii] <- FALSE
            
 
            
        }
        close(pb)
        ## merge back into data.frame
        chn <- cbind(terra::vect(gm[to_keep]),chn[to_keep,])
    }
    
    return(chn)
}



    ## if( any(chn$length < simplify_length) ){
    ##     message("Simplifying network")
    ##     gm <- terra::geom(chn) ## geom as a matrix to relabel
    ##     gnm <- colnames(gm)
    ##     ggm <- lapply(split(gm,gm[,'geom']),matrix,ncol=ncol(gm),col.names=gnm)
        
    ##     gm <- as.data.frame(gm) ## to make splitting easier
        
    ##     npart <- tapply(gm[,"part"],gm[,"geom"],max)
    ##     chn <- as.data.frame(chn) ## drop geom data

    ##     ## find channels to merge and order
    ##     idx <- which( chn$length < simplify_length )
    ##     idx <- idx[order(chn$length[idx])]

    ##     cnt <- 0
    ##     n <- length(idx)
    ##     print(paste("Evaluating",n,"potential candidates"))
    ##     for(ii in idx[1:10]){
    ##         if( chn$length[ii] >= simplify_length ){ next } ## catch incase a merge has already made it long
            
    ##         in_hn <-  which( chn$endNode == chn$startNode[ii] )
    ##         out_hn <-  which( chn$startNode == chn$startNode[ii] )
    ##         in_en <- which( chn$endNode == chn$endNode[ii] )
    ##         out_en <- which( chn$startNode == chn$endNode[ii] )
    ##         jj <- NA
    ##         if( length(in_en)==1 &
    ##             length(out_en)==1 ){
    ##             ## endNode is a 1-to-1 join and can merge d/s
    ##             jj <- out_en
    ##             chn$startNode[jj] <- chn$startNode[ii]
    ##         }
    ##         if( length(in_hn)==1 &
    ##             length(out_hn)==1 &
    ##             is.na(jj) ){
    ##             ## startNode is a 1-to-1 join and can merge u/s
    ##             jj <- in_hn
    ##             chn$endNode[jj] <- chn$endNode[ii]
    ##         }
    ##         if(is.na(jj)){ next } ## can't simplify
            
    ##         ## merge properties (default to those for jj)
    ##         chn$length[jj] <- chn$length[jj] + chn$length[ii]
    ##         chn$area[jj] <- chn$area[jj] + chn$area[ii]
    ##         chn$width[jj] <- chn$area[jj] / chn$length[jj]
    ##         chn$channelVol[jj] <- chn$channelVol[jj] + chn$channelVol[ii]
    ##         chn[ii,] <- NA

    ##         ## merge geom
    ##         kk <- npart[jj] + 1
    ##         kdx <-  gm[,"geom"]==ii
    ##         gm[kdx,"part"] <- kk
    ##         gm[kdx,"geom"] <- jj
    ##         npart[jj] <- kk
            
    ##         cnt <- cnt+1
    ##         #print(cnt)
    ##     }
        
    ##     ## merge back into data.frame
    ##     browser()
    ##     idx <- sort(unique(gm[,"geom"]))
    ##     gm[,"geom"] <- match(gm[,"geom"],idx)
    ##     chn <- terra::vect(gm,"polygons",atts=chn[idx,])
