#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
#' @examples
#' ## The vignettes contains more examples of the method calls.
#' 
#' ## create temporary directory for output
#' demo_dir <- tempfile("dygis")
#' dir.create(demo_dir)
#'
#' ## initialise processing
#' ctch <- dynatopGIS$new(file.path(demo_dir,"test"))
#'
#' ## add a catchment outline based on the digital elevation model
#' dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
#' dem <- terra::rast(dem_file)
#' dem <- terra::extend(dem,1)
#' catchment_outline <- terra::ifel(is.finite(dem),1,NA)
#' ctch$add_catchment(catchment_outline)
#' 
#' ## add digital elevation and channel data
#' ctch$add_dem(dem)
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' sp_lines <- terra::vect(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(sp_lines,property_names)
#' ctch$add_channel(chn)
#'
#' ## compute properties 
#' ctch$sink_fill() ## fill sinks in the catchment and computes dem flow directions
#' \donttest{
#' ctch$compute_properties() # like topograpihc index and contour length
#' ctch$compute_band()
#' ctch$compute_flow_lengths()
#' }
#' ## classify and create a model
#' \donttest{
#' ctch$classify("atb_20","atb",cuts=20) # classify using the topographic index
#' ctch$get_method("atb_20") ## see the details of the classification
#' ctch$combine_classes("atb_20_band",c("atb_20","band")) ## combine classes
#' ctch$create_model(file.path(demo_dir,"new_model"),"atb_20_band") ## create a model
#' list.files(demo_dir,pattern="new_model*") ## look at the output files for the model
#' }
#' ## tidy up
#' unlink(demo_dir)
#' @export
dynatopGIS <- R6::R6Class(
    "dynatopGIS",
    public = list(
        #' @description Initialise a project, or reopen an existing project
        #'
        #' @param projectFolder folder for data files
        #'
        #' @details This loads the project data files found in the \code{projectFolder} if present. If not the folder is created. The project data files are given by \code{projectFolder/<filename>.<tif,shp>}
        #'
        #' @return A new `dynatopGIS` object
        initialize = function(projectFolder){
            private$apply_initialize( projectFolder )
            invisible(self)
        },                
        #' @description Add a catchment outline to the `dynatopGIS` project
        #'
        #' @param catchment a \code{SpatRaster} object or the path to file containing one which contains a rasterised catchment map.
        #'
        #' @details If not a \code{SpatRaster} object the the catchment is read in using the terra package. Finite values in the raster indicate that the area is part of the catchment; with each subcatchment taking a unique finite value. Note that in the later processing it is assumed that outflow from the subcatchments can occur only through the channel network. The resolution and projection of the project is taken from the provided catchment
        #'
        #' @return \code{invisible(self)}
        add_catchment = function(catchment){
            if(!("SpatRaster" %in% class(catchment))){ catchment <- terra::rast(as.character(catchment)) }
            if(!("SpatRaster" %in% class(catchment))){ stop("catchment is not a SpatRaster") }
            private$apply_add_catchment(catchment)
            invisible(self)
        },
        #' @description Import a dem to the `dynatopGIS` object
        #'
        #' @param dem a \code{raster} layer object or the path to file containing one which is the DEM
        #' @param fill_na  should NA values in dem be filled. See details
        #' @param verbose Should additional progress information be printed
        #'
        #' @details If not a \code{raster} the DEM is read in using the terra package. If \code{fill_na} is \code{TRUE} all NA values other then those that link to the edge of the dem are filled so they can be identified as sinks.
        #'
        #' @return suitable for chaining              
        add_dem = function(dem,fill_na=-9999){
            if(!("SpatRaster" %in% class(dem))){ dem <- terra::rast(as.character(dem)) }
            if(!("SpatRaster" %in% class(dem))){ stop("dem is not a SpatRaster") }
            private$apply_add_dem(dem,fill_na)
            invisible(self)
        },
        #' @description Import channel data to the `dynatopGIS` object
        #'
        #' @param channel a SpatVect object or file path that can be loaded as one containing the channel information
        #' @param verbose Should additional progress information be printed
        #' @details Takes the representation of the channel network as a SpatVect with properties name, length, area, startNode, endNode and overlaying it on the DEM. In doing this a variable called id is created (or overwritten) other variables in the data frame are passed through unaltered.
        #'
        #' @return suitable for chaining
        add_channel = function(channel,verbose=FALSE){
            if(!is(channel,"SpatVector")){ channel <- terra::vect( as.character(channel) ) }
            if(!is(channel,"SpatVector")){ stop("channel is not a SpatVector object") }

            private$apply_add_channel(channel,as.logical(verbose))
            invisible(self)
        },        
        #' @description Add a layer of geographical information
        #'
        #' @param layer the raster layer to add (see details)
        #' @param layer_name name to give to the layer
        #'
        #' @details The layer should either be a raster layer or a file that can be read by the \code{raster} package. The projection, resolution and extent are checked against the existing project data. Only layer names not already in use (or reserved) are allowed. If successful the layer is added to the project tif file.
        #' @return suitable for chaining
        add_layer = function(layer,layer_name=names(layer)){
            layer_name <- as.character(layer_name)
            if(!("SpatRaster" %in% class(layer))){ layer <- terra::rast(as.character(layer)) }
            if(!("SpatRaster" %in% class(layer))){ stop("layer is not a SpatRaster") }
            if( length(layer_name) != terra::nlyr(layer) ){ stop("Length of names does not match number of layers") }
            private$apply_add_layer(layer,layer_name)
            invisible(self)
        },
        #' @description Get a layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a `raster` layer of the requested information if layer_name is given else a vector of layer names
        get_layer = function(layer_name=character(0)){
            ## create a data frame of available layers
            tmp <- data.frame(layer = names(private$brk), source = terra::sources(private$brk))
            if( "channel" %in% tmp$layer ){ tmp <- rbind(tmp, data.frame(layer="channel_vect", source= terra::sources(private$shp))) }
            
            ## handle case where a list of layers is requested
            if( length(layer_name) == 0 ){ return(tmp) }

            ## check layer name exists
            layer_name <- match.arg(layer_name,tmp$layer)
            
            ## make raster and return
            if( layer_name == "channel_vect" ){
                return( private$shp )
            }else{
                return( private$brk[[layer_name]] )
            }
        },
        #' @description Plot a layer
        #' @param layer_name the name of layer to plot
        #' @param add_channel should the channel be added to the plot
        #' @return a plot
        plot_layer = function(layer_name,add_channel=TRUE){
            lyr <- self$get_layer(layer_name)
            terra::plot( lyr, main = layer_name)
            if( add_channel & length(private$shp) > 0){
                terra::plot(private$shp, add=TRUE )
            }
        },
        #' @description The sink filling algorithm of Planchona and Darboux (2001)
        #'
        #' @param min_grad Minimum gradient between cell centres
        #' @param max_it maximum number of replacement cycles
        #' @param verbose print out additional diagnostic information
        #' @param hot_start start from filled_dem if it exists
        #' @param flow_type The type of flow routing to apply see details
        #' @details The algorithm implemented is based on that described in Planchona and Darboux, "A fast, simple and versatile algorithm to fill the depressions in digital elevation models" Catena 46 (2001). A pdf can be found at (<https://horizon.documentation.ird.fr/exl-doc/pleins_textes/pleins_textes_7/sous_copyright/010031925.pdf>). The adaptations made are to ensure that all cells drain only within the subcatchments if provided.
        #'
        #' The flow_type can be either
        #' - "quinn" where flow is split across all downslope directions or
        #' - "d8" where all flow follows the steepest between cell gradient
        #'
        sink_fill = function(min_grad = 1e-4,max_it=1e6,verbose=FALSE, hot_start=FALSE, flow_type=c("quinn","d8")){
            flow_type <- match.arg(flow_type)
            private$apply_sink_fill(min_grad,max_it,verbose,hot_start,flow_type)
            invisible(self)
        },
        #' @description Computes the computational band of each cell
        #'
        #' @param type type of banding
        #' @param verbose print out additional diagnostic information
        #'
        #' @details Banding is used within the model to define the HRUs and control the order of the flow between them; HRUs can only pass flow to HRUs in a lower numbered band. Currently only a strict ordering of river channels and cells in the DEM is implemented. To compute this the algorithm passes first up the channel network (with outlets being in band 1) then through the cells of the DEM in increasing height.
        compute_band = function(type=c("strict"), verbose=FALSE){
            type = match.arg(type)
            private$apply_band(type,verbose)
            invisible(self)
        },
        #' @description Computes statistics e.g. gradient, log(upslope area / gradient) for raster cells
        #'
        #' @param min_grad gradient that can be assigned to a pixel if it can't be computed
        #' @param verbose print out additional diagnostic information
        #'
        #' @details The algorithm passed through the cells in decreasing height. Min grad is applied to all cells. It is also used for missing gradients in pixels which are partially channel but have no upslope neighbours.
        compute_properties = function(min_grad = 1e-4,verbose=FALSE){
            private$apply_compute_properties(min_grad,verbose)
            invisible(self)
        },
        #' @description Computes flow length for each pixel to the channel
        #'
        #' @param flow_routing TODO
        #' @param verbose print out additional diagnostic information
        #'
        #' @details The algorithm passes through the cells in the DEM in increasing height. Three measures of flow length to the channel are computed. The shortest length (minimum length to channel through any flow path), the dominant length (the length taking the flow direction with the highest fraction for each pixel on the path) and expected flow length (flow length based on sum of downslope flow lengths based on fraction of flow to each cell). By definition cells in the channel that have no land area have a length of NA.
        compute_flow_lengths = function(flow_routing=c("expected","dominant","shortest"), verbose=FALSE){
            flow_routing = match.arg(flow_routing)
            private$apply_flow_lengths(flow_routing,verbose)
            invisible(self)
        }, 
        #' @description Create a catchment classification based cutting an existing layer into classes
        #' @param layer_name name of the new layer to create
        #' @param base_layer name of the layer to be cut into classes
        #' @param cuts values on which to cut into classes. These should be numeric and define either the number of bands (single value) or breaks between band (multiple values).
        #'
        #' @details This applies the given cuts to the supplied landscape layer to produce areal groupings of the catchment. Cuts are implement using \code{terra::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        classify = function(layer_name,base_layer,cuts){
            private$apply_classify(as.character(layer_name), as.character(base_layer), cuts)
            invisible(self)
        },
        #' @description Combine any number of classifications based on unique combinations and burns
        #' @param layer_name name of the new layer to create
        #' @param pairs a vector of layer names to combine into new classes through unique combinations. Names should correspond to raster layers in the project directory.
        #' @param burns a vector of layer names which are to be burnt on
        #'
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given. Cuts are implement using \code{terra::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        combine_classes = function(layer_name,pairs,burns=NULL){
            private$apply_combine_classes(as.character(layer_name),
                                          as.character(pairs),
                                          as.character(burns))
            invisible(self)
        },
        #' @description Compute a Dynamic TOPMODEL
        #'
        #' @param layer_name name for the new model and layers
        #' @param class_layer the layer defining the topographic classes
        #' @param sf_opt Surface solution to use
        #' @param sz_opt transmissivity profile to use
        #' @param rain_layer the layer defining the rainfall inputs
        #' @param rain_label Prepended to rain_layer values to give rainfall series name        
        #' @param pet_layer the layer defining the pet inputs
        #' @param pet_label Prepended to pet_layer values to give pet series name
        #' @param verbose print more details of progress
        #'
        #' @details The \code{class_layer} is used to define the HRUs. Flow between HRUs is based on the ordering of the catchment (see the \code{compute_band} method). Flow from a HRU can only go to a HRU with a lower band.
        #' Setting the sf_opt and sz_opt options ensures the model is set up with the correct parameters present.
        #' The \code{rain_layer} (\code{pet_layer}) can contain the numeric id values of different rainfall (pet) series. If the value of \code{rain_layer} (\code{pet_layer}) is not \code{NULL} the weights used to compute an averaged input value for each HRU are computed, otherwise an input table for the models generated with the value "missing" used in place of the series name.
        create_model = function(layer_name,class_layer,
                                sf_opt = c("cnst","kin"),
                                sz_opt = c("exp","bexp","cnst","dexp"),
                                rain_layer=NULL, rain_label=character(0),
                                pet_layer=NULL, pet_label=character(0),
                                verbose=FALSE){
            
            ## check valid transmissivity and channel_solver
            sf_opt<- match.arg(sf_opt)
            sz_opt <- match.arg(sz_opt)
            
            private$apply_create_model(class_layer,
                                       rain_layer, rain_label,
                                       pet_layer, pet_label,
                                       layer_name,verbose,
                                       sf_opt, sz_opt)

            invisible(self)
        },
        #' @description locate gauges on the channel newtork
        #' @param gauge_map a SpatVect defining the gauge locations
        #' @param gauge_name feild in SpatVect containing the name of the gauge
        locate_gauges = function(gauge_map,gauge_name="name"){
            if(!("SpatVect" %in% class(gauge_map))){ gauge_map <- terra::vect(as.character(gauge_map)) }
            if(!("Spatvect" %in% class(gauge_map))){ stop("gauge_map is not a SpatVect") }
            if(!(gauge_name %in% names(gauge_map))){ stop("gauge_name is not a proerty of gauge_map") }

            private$apply_locate_gauges(gauge_map,gauge_name)
        },
        #' @description get the version number
        #' @return a numeric version number
        #' @details the version number indicates the version of the algorithms within the object
        get_version = function(){
            private$version
        },
        #' @description get the cuts and burns used to classify
        #' @param layer_name the name of layer whose classification method is returned
        #' @return a list with two elements, cuts and burns
        get_method = function(layer_name){
            ## check layer name exists
            layer_name <- match.arg(layer_name,names(private$brk))
            
            jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[layer_name]])),".json")
            if( !file.exists(jsonFile) ){
                stop("No json file giving basis of the classifications")
            }

            return( jsonlite::fromJSON( jsonFile ) )
        }               
    ),
    private = list(
        version = "0.3.1",
        projectFolder = character(0),
        brk = character(0),
        shp = character(0),
        reserved_layers = c("catchment","dem","channel","filled_dem",
                            "gradient","upslope_area","atb",
                            "band","flow_length"),
        readJSON = function(fn){
            
            if(!file.exists(fn)){
                return(NULL) #warning("No method json file found")
            }
            jsonlite::fromJSON(fn)
        },
        ## check and read the project files
        apply_initialize = function(projectFolder){
            ## see if folder exists
            if( !dir.exists(projectFolder) ){
                message( "Creating new folder" )
                dir.create(projectFolder, showWarnings = TRUE, recursive = TRUE)
            }

            ## read in raster data
            rstNames <- list.files(projectFolder, pattern=".tif$",full.names=TRUE)
            if(length(rstNames) == 0){
                brk <- chn <- NULL
                message( paste("Starting new project at", projectFolder) )
            }else{
                message( paste("Reading existing project at", projectFolder) )
                brk <- terra::rast(list.files(projectFolder, pattern=".tif$",full.names=TRUE))

                if( !all.equal(diff(terra::res(brk)),0) | terra::is.lonlat(brk) ){
                    stop("Raster data not valid: Processing currently only works on projected data with a square grid")
                }
                if( !("catchment" %in% names(brk)) ){ stop("No catchment layer in Raster data") }
                ## read in shape file and check
                if( file.exists( file.path(projectFolder,"channel.shp")) ){
                    chn <- terra::vect( file.path(projectFolder,"channel.shp") )
                
                    if(  terra::crs(brk, proj=TRUE) != terra::crs(chn, proj=TRUE) ){
                        stop("Shape and Raster projections do not match")
                    }
                    if( !("channel" %in% names(brk)) ){
                        stop("No channel layer in the raster file but channel shapefile specified")
                    }
                    if( !file.exists( file.path(projectFolder,"channel.rds") ) ){
                        stop("No channel flow direction file but channel shapefile specified")
                    }
                    
                }else{
                    chn <- NULL
                    if( ("channel" %in% names(brk)) ){
                        stop("Channel layer in the raster file but channel shapefile missing")
                    }
                }
            }

            
            
            private$projectFolder <- projectFolder
            private$brk <- brk
            private$shp <- chn
        },
        ## add a catchment map
        apply_add_catchment = function(catchment){
            
            if("catchment" %in% names(private$brk) ){
                stop("The catchment map already exists, start a new project")
            }
            
            ## check the projection is OK
            if( !all.equal(diff(terra::res(catchment)),0) | terra::is.lonlat(catchment) ){
                stop("Processing currently only works on projected maps with a square grid")
            }

            ## check it is na padded
            if( any( is.finite(catchment[1,][[1]]) ) |
                any( is.finite(catchment[nrow(catchment),][[1]]) ) |
                any( is.finite(catchment[,1][[1]]) ) |
                any( is.finite(catchment[,ncol(catchment)][[1]]) ) ){
                stop("catchment must be padded with NA rows and columns")
            }

            ## save
            names(catchment) <- "catchment"
            demFile <- file.path(private$projectFolder,"catchment.tif")
            terra::writeRaster(catchment, demFile)
            private$brk <- terra::rast(demFile)
        },
        ## adding dem
        apply_add_dem = function(dem,fill_na){

            if(!("catchment" %in% names(private$brk)) ){
                stop("The catchment map does not exists, try running add_catchment")
            }
            
            if("dem" %in% names(private$brk) ){
                stop("The DEM already exists")
            }

            if( !terra::compareGeom(dem,private$brk,stopOnError=FALSE) ){
                stop("New layer does not match resolution, extent or projection of project")
            }

            ## convert na values
            dem <- terra::ifel(is.na(dem), fill_na, dem)
            
            ## mask layer to catchment
            dem <- terra::mask(dem,private$brk[[ "catchment" ]])

            ## save
            names(dem) <- "dem"
            rstFile <- file.path(private$projectFolder,"dem.tif")
            terra::writeRaster(dem, rstFile)
            private$brk <- c(private$brk, terra::rast( rstFile ))

        },
        ## add the channel
        apply_add_channel = function(chn,verbose){

            rq <- c("catchment")
            if(!all(rq %in% names(private$brk))){
                stop("Not all required layers are available")
            }
            
            if( length(private$shp) > 0 ){
                stop("Channel already added")
            }
            
            ## check the required field names are present
            if( !all( c("name","length","area","startNode","endNode","width","slope") %in% names(chn)) ){
                stop("A required property name is not specified")
            }
            
            ## check projection of the chn
            if( terra::crs(private$brk, proj=TRUE) != terra::crs(chn, proj=TRUE) ){
                stop("Projection of channel object does not match that of project")
            }
            
            ## check if there is an id feild which will be overwritten
            if( ("id" %in% names(chn)) ){
                warning("The name id is reserved and will be overwritten in the channel import")
            }

            ## ensure required properties are of correct type
            chn$name <- as.character(chn$name)
            chn$length <- as.numeric(chn$length)
            chn$area <- as.numeric(chn$area)
            chn$startNode <- as.character(chn$startNode)
            chn$endNode <- as.character(chn$endNode)
            chn$width <- as.numeric(chn$width)
            chn$slope <- as.numeric(chn$slope)
            
            if( !all(is.finite(chn$length)) ){ stop("Some non-finite values of length found!") }
            if( !all(is.finite(chn$width)) ){ stop("Some non-finite values of width found!") }
            if( !all(is.finite(chn$slope)) ){ stop("Some non-finite values of slope found!") }

            ## TODO - the whole of the following can probably be simplified....


            ## arrange id and band in order of flow direction - so lowest values at outlets of the network
            if( verbose ){ print("Computing channel id's and bands") }
            ## This is much quicker using vectors and not constantly accessing via the vect object
            id <- rep(as.integer(0),nrow(chn))
            bnd <- rep(as.integer(0),nrow(chn))
            sN <- chn$startNode
            eN <- chn$endNode
          
            idx <- !(eN %in%sN) ## outlets are channel lengths whose outlet does not join another channel
            it <- 1
            cnt <- table(sN) ## we should never vist a node more times then it is a starting point
            while(sum(idx)>0){
                id[idx] <- max(id) + 1:sum(idx)
                bnd[idx] <- it
                
                jdx <- sN[idx]
                for(ii in jdx){ cnt[ii] <- cnt[ii] - 1 } ## since sN might appear more then once..
                if( any(cnt<0) ){
                    stop(paste("Loop involving nodes:", paste(names(cnt)[cnt<0],collapse=", ")))
                }
                jdx <- jdx[cnt[jdx]==0] ## only move up if it is the last visit to the startNode
                idx <- eN %in% jdx
                it <- it+1
            }
            chn$id <- id
            chn$band <- bnd
            if( !(all(chn$id > 0))){ stop("Error ingesting channel: check connectivity") }
            if( !(all(chn$band > 0))){ stop("error ingesting channel: problem with bands") }
            if( !(all(cnt==0)) ){ stop("error ingesting channel: problem with visiting all points") }
            chn <- chn[ order(chn$id),]

            ## create a raster of channel id numbers
            ## TODO - possibly sort so do biggest area first???
            chn_rst <- terra::rasterize(chn,private$brk[["catchment"]],field = "id",touches=TRUE)
            chn_rst <- terra::mask(chn_rst,private$brk[["catchment"]]) ## make sure value occur only on cells within the catchment
            names(chn_rst) <- "channel"

            ## compute channel routing - order flow links from bottom to top
            if( verbose ){ print("Computing channel routing") }
            chn_route <- rep(list(NULL),nrow(chn))
            id <- chn$id
            sN <- chn$startNode
            eN <- chn$endNode
            for(ii in 1:nrow(chn)){
                idx <- which( sN==eN[ii] )
                if( length(idx) == 0 ){ next }
                chn_route[[ii]] <- cbind( id[ii], id[idx], 1/length(idx) )
            }
            chn_route <- do.call(rbind,chn_route)
            colnames(chn_route) <- c("from","to","fraction")
            if( !all( chn_route[,"from"] > chn_route[,"to"] ) ){ stop("Incorrect channel routing") }

            ## save output
            rdsFile <- file.path(private$projectFolder,"channel.rds")
            shpFile <- file.path(private$projectFolder,"channel.shp")
            rstFile <- file.path(private$projectFolder,"channel.tif")
            saveRDS(chn_route,rdsFile)
            terra::writeVector(chn, shpFile)
            terra::writeRaster(chn_rst,rstFile, names = "channel")
            
            private$brk <- c(private$brk, terra::rast(rstFile))
            private$shp <- terra::vect(shpFile)
        },
        ## Add a layer
        apply_add_layer=function(layer,layer_name){          
            if( any(layer_name %in% private$reserved_layers) ){
                stop("Name is reserved")
            }
            if( any(layer_name %in% names(private$brk)) ){
                stop("Name is already used")
            }
            if( !terra::compareGeom(layer,private$brk,stopOnError=FALSE) ){
                ## try buffering it as for dem when read in
                layer <- terra::extend(layer,c(1,1))
            }
            if( !terra::compareGeom(layer,private$brk,stopOnError=FALSE) ){
                stop("New layer does not match resolution, extent or projection of project")
            }
            names(layer) <- layer_name
            rstFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            terra::writeRaster(layer, rstFile)
            private$brk <- c(private$brk, terra::rast( rstFile ))
         },
        ## Sink fill
        apply_sink_fill = function(min_grad,max_it,verbose,hot_start,flow_type){
            ## recall the catchments is padded with NA values
            
            rq <- ifelse(hot_start,
                         c("filled_dem","channel","catchment"),
                         c("dem","channel","catchment"))
            if(!all(rq %in% names(private$brk))){
                stop("Not all required layers are available")
            }

            d <- ifelse(hot_start,"filled_dem","dem")
            d <- terra::as.matrix( private$brk[[d]] ,wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE )
            ctch <- terra::as.matrix( private$brk[["catchment"]] , wide=TRUE )
            
            ## values that should be valid
            to_be_valid <- !is.na(ctch) #!is.na(d) & is.na(ch)  # all values not NA should have a valid height
            is_valid <- !is.na(ch) # & !is.na(d) # TRUE if a channel cell for initialisation
            changed <- is_valid # cells changed at last iteration
            fd <- to_be_valid*Inf; fd[is_valid] <- d[is_valid]
            to_eval <- is_valid; to_eval[] <- FALSE

            ## dimensions
            nf <- sum(!is.na(ctch))
            
            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- matrix(sqrt(sum(rs^2)),3,3)
            dxy[1,2] <- dxy[3,2] <- rs[2]
            dxy[2,1] <- dxy[2,3] <- rs[1]
            dxy[2,2] <- 0
            dxy <- min_grad*dxy
            
            it <- 1
            ## start of iteration loop
            while(any(changed[]) & it<=max_it){

                ## work out all the the cells to evaluate
                ## should be next to those that are changed
                to_eval[] <- FALSE
                idx <- which(changed,arr.ind=TRUE) # index of changed cells
                sctch <- ctch[idx] # sub catchment number
                ##if( any(is.na(sctch)) ){ browser() }
                jdx <- idx
                for(ii in c(-1,0,1)){
                    for(jj in c(-1,0,1)){
                        if(ii==0 & jj==0){next}
                        ## adjust jdx
                        ##browser()
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        ## trim so only evaluate
                        ## (ii) cells within the same subcatchment
                        cjdx <- ctch[jdx]
                        kdx <- !is.na(cjdx) & (sctch == cjdx)
                        jjdx <- jdx[kdx,, drop=F]
                        to_eval[jjdx] <- !is_valid[jjdx] #TRUE
                    }
                }
                to_eval <- to_eval & to_be_valid #& !is_valid
                if(verbose){
                    cat("Iteration",it,"\n")
                    cat("\t","Cells to evaluate:",sum(to_eval),"\n")
                    cat("\t","Percentage Complete:",
                        round(100*sum(is_valid)/nf,1),"\n") #to_be_valid),1),"\n")
                }
                ## alter min value for the evaluation cells
                idx <- which(to_eval,arr.ind=TRUE) # index of changed cells
                sctch <- ctch[idx] # sub catchment number
                ##if( any(is.na(sctch)) ){ browser() }
                jdx <- idx
                mind <- rep(Inf,nrow(idx))
                for(ii in c(-1,0,1)){
                    for(jj in c(-1,0,1)){
                        if(ii==0 & jj==0){next}
                        ## adjust jdx
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        ## trim so only evaluate
                        ## (ii) cells within the same subcatchment
                        cjdx <- ctch[jdx]
                        kdx <- !is.na(cjdx) & (sctch == cjdx)
                        if( !any(kdx) ){ next }
                        jjdx <- jdx[kdx,,drop=F]
                        ##if( any(is.na(kdx)) | any(is.na(jjdx)) ){ browser() }
                        ##if( sum(kdx)==0 ){ browser() }##| sum(kdx) != nrow(jjdx) ){ browser() }
                        tmp <- pmin(mind[kdx],fd[jjdx] + dxy[ii+2,jj+2],na.rm=TRUE)
                        ##if( length(tmp) != length( mind[kdx] )){ browser() }
                        mind[kdx] <- pmin(mind[kdx],fd[jjdx] + dxy[ii+2,jj+2],na.rm=TRUE)
                    }
                }
                changed[] <- FALSE
                is_valid[idx] <- d[idx]>mind ## cells where the dem value is valid
                mind <- pmax(d[idx],mind) ## mind is now the replacemnt value
                changed[idx] <- mind < fd[idx]
                fd[idx] <- mind
                
                ## end of loop
                it <- it+1
            }

            rfd <- terra::rast( private$brk[["dem"]], names="filled_dem", vals=fd )
            rstFile <- file.path(private$projectFolder,"filled_dem.tif")
            if(hot_start){
                terra::writeRaster(rfd, rstFile,overwrite=TRUE)
                private$brk[["filled_dem"]] <- terra::rast(rstFile)
            }else{
                terra::writeRaster(rfd, rstFile)
                private$brk <- c( private$brk, terra::rast(rstFile))
            }
     
            if(it>max_it){ stop("Maximum number of iterations reached, sink filling not complete") }
            
            ## work out routing
            ## private$apply_flow_direction(flow_type,verbose)
        },
        ## Function to compute the bands
        apply_band = function(type,verbose){
            rq <- c("filled_dem","channel","catchment")
            if(!all( rq %in% names( private$brk) )){
                stop("Not all required input layers have been generated \n",
                     "Try running sink_fill first")
            }
            
            ## rasterize channel band to start
            rbnd <- terra::rasterize(private$shp, private$brk[["catchment"]],field = "band",touches=TRUE)
            
            ## load raster layer
            d <- terra::as.matrix( private$brk[["filled_dem"]], wide=TRUE )
            bnd <- terra::as.matrix( rbnd,  wide=TRUE )
            ctch <- terra::as.matrix( private$brk[["catchment"]],  wide=TRUE )
            
            if( verbose ){ print("Computing hillslope routing and band") }

            ## distances and contour lengths
            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## if we go up in height order then we are working from near the channel to the heighest point
            idx <- order(d,na.last=NA)
            
            ## create flow direction storage
            fd <- matrix(as.numeric(NA),length(idx),11)
            colnames(fd) <- c("cell","gcl","dcl","topLeft","left","bottomLeft","top","bottom","topRight","right","bottomRight")
            fd_cnt <- 0
            
            n_to_eval <- length(idx)
            
            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }else{
                next_print <- Inf
            }
            
            w <- rep(0,8)
            for(ii in idx){
                ## skip if in a channel
                if(is.finite( bnd[ii])){ next }
                
                ## it is not a channel
                w[] <- 0
                jdx <- ii+delta
                grd <- (d[jdx]-d[ii])/dxy
                gcl <- grd*dcl
                cjdx <- ctch[jdx]
                is_lower <- is.finite(gcl) & gcl<0 & is.finite(cjdx) & cjdx==ctch[ii]
                ## compute weights
                if(type == "d8"){
                    grd[!is_lower] <- Inf
                    is_lower[] <- FALSE
                    is_lower[ which.min(grd) ] <- TRUE
                }
                sum_gcl <- sum( gcl[is_lower] )
                sum_dcl <- sum( dcl[is_lower] )
                w[is_lower] <- gcl[is_lower] / sum_gcl
                
                fd_cnt <- fd_cnt + 1
                fd[fd_cnt,] <- c(ii,sum_gcl,sum_dcl,w)

                ## compute the band
                bnd[ii] <- max( bnd[jdx[w>0]] ) + 1
                
                ## verbose output here
                if(it >= next_print){
                    cat(round(100*it / n_to_eval,1),
                        "% complete","\n")
                    next_print <- next_print+print_step
                }
                
                it <- it+1
            }

            ## save flow direction
            saveRDS(fd[1:fd_cnt,],file.path(private$projectFolder,"dem.rds"))

            ## save band
            terra::values(rbnd) <- bnd
            names(rbnd) <- "band"
            rstFile <- file.path(private$projectFolder,"band.tif")
            terra::writeRaster(rbnd, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))
        },
        
        ## ## work out banding of the cells
        ## apply_band = function(type,verbose){
            
        ##     rq <- c("filled_dem","channel")
        ##     if(!all( rq %in% names( private$brk) )){
        ##         stop("Not all required input layers have been generated \n",
        ##              "Try running sink_fill first")
        ##     }
            
        ##     ## load raster layer
        ##     ch <- terra::as.matrix( private$brk[["channel"]],  wide=TRUE )
            
        ##     rq <- c( file.path(private$projectFolder,"dem.rds"),
        ##             file.path(private$projectFolder,"channel.rds") )
        ##     if( ! all( file.exists(rq) ) ){
        ##         stop("No flow routing records defined\n",
        ##              "Try running compute_flow_paths first")
        ##     }else{
        ##         flow_routing <- readRDS(rq[1])
        ##         channel_routing <- readRDS( rq[2] )
        ##     }
            
        ##     ## compute the channel band
        ##     if( verbose ){ print("Computing channel bands") }
        ##     id <- private$shp$id
        ##     bnd <- rep(-999,length(id))
        ##     idx <- !(id %in% channel_routing[,"from"]) ## locations with no flow from them
        ##     bnd[idx] <- 1
        ##     for(ii in 1:nrow(channel_routing)){
        ##         bnd[channel_routing[ii,"from"]] <- max( bnd[channel_routing[ii,"from"]],
        ##                                                bnd[channel_routing[ii,"to"]] + 1 )
        ##     }
            
            
        ##     ## convert ch to bnd matrix
        ##     if(verbose){ print("Making band matrix") }
        ##     bnd <- bnd[match(ch,id)]
        ##     dim(bnd) <- dim(ch)
            
            
        ##     ## set up loop through hillslope
        ##     if(verbose){ print("Setting up hillslope loop") }
            
        ##     n_to_eval <- nrow(flow_routing)
            
        ##     nr <- nrow(ch); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
        ##     it <- 1
        ##     if(verbose){
        ##         print_step <- round(n_to_eval/20)
        ##         next_print <- print_step
        ##     }else{
        ##         next_print <- Inf
        ##     }
            
        ##     ## loop from the lowest link
        ##     if(verbose){ print("Starting hillslope loop") }
        ##     for(rw in 1:nrow(flow_routing)){
        ##         ii <- flow_routing[rw,1]
        ##         w <- flow_routing[rw,-1]
        ##         jj <- ii + delta
        ##         jj <- jj[w>0]
        ##         bnd[ii] <- max( bnd[jj] )+1
                
        ##         ## verbose output here
        ##         if(it >= next_print){
        ##             cat(round(100*it / n_to_eval,1),
        ##                 "% complete","\n")
        ##             next_print <- next_print+print_step   
        ##         }
                
        ##         it <- it+1
        ##     }
            
        ##     ## save raster maps
        ##     out <- terra::rast( private$brk[["dem"]], names="band", vals=bnd)
        ##     rstFile <- file.path(private$projectFolder,"band.tif")
        ##     terra::writeRaster(out, rstFile); 
        ##     private$brk <- c( private$brk, terra::rast(rstFile))
            
        ##     shpFile <- file.path(private$projectFolder,"channel.shp")
        ##     terra::writeVector(private$shp, shpFile, overwrite=TRUE)
        ## },
        
        ## Function to compute the properties
        apply_compute_properties = function(min_grad,verbose){
           
            if( verbose ){ print("Loading data") }
            rq <- c("filled_dem","channel")
            if(!all( rq %in% names( private$brk) )){
                stop("Not all required input layers have been generated \n",
                     "Try running sink_fill first")
            }

            rq <- c( file.path(private$projectFolder,"dem.rds"),
                    file.path(private$projectFolder,"channel.rds") )
            if( ! all( file.exists(rq) ) ){
                stop("No flow routing records defined\n",
                     "Try running compute_flow_paths first")
            }else{
                flow_routing <- readRDS(rq[1])
                channel_routing <- readRDS( rq[2] )
            }
            
            ## load raster layer
            d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE)
            
            if( verbose ){ print("Setting up output") }

            ## work out order to pass through the cells
            n_to_eval <- nrow(flow_routing)
            
            ## distance between cell centres
            rs <- terra::res( private$brk )
            ##dxy <- rep(sqrt(sum(rs^2)),8)
            ##dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            ##dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## initialise output
            gr <- upa <- atb <- d*NA
            upa[is.finite(d)] <- prod(rs) ## initialise upslope area from resolution
            
            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }else{
                next_print <- Inf
            }

            if( verbose ){ print("Computing hillslope") }
            ## loop downslope
            for(rwnum in nrow(flow_routing):1){
                ii <- flow_routing[rwnum,1]
                gcl <- flow_routing[rwnum,2] # sum of gradient * contour length for all d/s
                dcl <- flow_routing[rwnum,3] # sum of contour length for all d/s
                w <- flow_routing[rwnum,4:11] # weight of flow direction
                to_use <- w>0
                
                if( !is.na(ch[ii]) ){ stop("Propergating from a channel cell") }
                
                if( !any(w>0) ){ stop(paste("Cell",ii,"is a hillslope cell with no outflows")) }
                
                ngh <- ii + delta ## neighbouring cells
                ##grd <- (d[ii]-d[ngh])/dxy ## compute gradient

                ##gcl <- grd[to_use]*dcl[to_use]
                ## gradient
                ##gr[ii] <- max(sum(gcl) / sum(dcl[to_use]),min_grad)
                gr[ii] <- max( gcl/dcl, min_grad )
                ## topographic index
                atb[ii] <- log(upa[ii]/gr[ii]) #log( upa[ii] / sum(gcl) )
                ## propogate area downslope
                upa[ ngh ]  <- upa[ ngh ] + w*upa[ii]
                                
                ## verbose output here
                if(it >= next_print){
                    cat(round(100*it / n_to_eval,1),
                        "% complete","\n")
                    next_print <- next_print+print_step
                }
                
                it <- it+1
            }

            ## merge upslope areas into the channel object
            ch_upa <- tapply(upa,ch,sum)
            ch_upa <- ch_upa[setdiff(names(ch_upa),"NaN")]
            idx <- match(paste(private$shp$id),names(ch_upa))
            private$shp$up_area = as.numeric(NA)
            private$shp$up_area[idx] <- as.numeric(ch_upa)
            if( !all(is.finite(private$shp$up_area)) ){ stop("All upslope channel areas should be finite") }

            ## remove channel area bit from upa
            upa <- upa * !is.finite(ch)
            
            ## compute catchment area to each reach
            ct_area <- private$shp$up_area
            for(ii in nrow(channel_routing):1){
                ct_area[channel_routing[ii,"to"]] <- ct_area[channel_routing[ii,"to"]] +
                    ct_area[channel_routing[ii,"from"]]*channel_routing[ii,"fraction"]
            }
            private$shp$ct_area <- ct_area
            
            ## save raster maps
            out <- terra::rast( private$brk[["dem"]], names="gradient", vals=gr )
            rstFile <- file.path(private$projectFolder,"gradient.tif")
            terra::writeRaster(out, rstFile);
            private$brk <- c( private$brk, terra::rast(rstFile))

            out <- terra::rast( private$brk[["dem"]], names="upslope_area", vals=upa )
            rstFile <- file.path(private$projectFolder,"upslope_area.tif")
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))

            out <- terra::rast( private$brk[["dem"]], names="atb", vals=atb )
            rstFile <- file.path(private$projectFolder,"atb.tif")
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))

            shpFile <- file.path(private$projectFolder,"channel.shp")
            terra::writeVector(private$shp, shpFile, overwrite=TRUE)

        },
        ## work out flow lengths to channel
        apply_flow_lengths = function(type,verbose){

            type <- paste0(type,"_flow_length")
            
            ## check not already computed
            if(type %in% names(private$brk)){ stop("Already computed") }
            
            rq <- c("channel")
            if(!all( rq %in% names( private$brk) )){
                stop("Not all required input layers have been generated \n",
                     "Try running sink_fill first")
            }

            ## load raster layer
            ##d <- terra::as.matrix( private$brk[["filled_dem"]], wide=TRUE )
            ch <- terra::as.matrix( private$brk[["channel"]],  wide=TRUE )

            rq <- c( file.path(private$projectFolder,"dem.rds"),
                    file.path(private$projectFolder,"channel.rds") )
            if( ! all( file.exists(rq) ) ){
                stop("No flow routing records defined\n",
                     "Try running compute_flow_paths first")
            }else{
                flow_routing <- readRDS(rq[1])
                channel_routing <- readRDS( rq[2] )
            }
            
            ## create a distance matrix, initialise with channel elements 0
            fl <- ch; fl[fl>0] <- 0
            
            ## distances and contour lengths
            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            #dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(fl); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

            n_to_eval <- nrow(flow_routing)

            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }else{
                next_print <- Inf
            }

            w <- rep(0,8)
            tmp <- rep(NA,8)
            for(rw in 1:nrow(flow_routing)){
                ii <- flow_routing[rw,1]
                w[] <- flow_routing[rw,4:11]
                jj <- (ii + delta)
                tmp[] <- (fl[jj] + dxy)

                if(type=="shortest_flow_length"){
                    fl[ii] <- min( tmp[w>0] )
                }
                if(type=="dominant_flow_length"){
                    fl[ii] <- tmp[ which.max(w) ]
                }
                if(type=="expected_flow_length"){
                    fl[ii] <- sum( tmp[w>0] * w[w>0] )
                }
                
                ## verbose output here
                if(it >= next_print){
                    cat(round(100*it / n_to_eval,1),
                        "% complete","\n")
                    next_print <- next_print+print_step
                }
                
                it <- it+1
            }


            out <- terra::rast( private$brk[["dem"]], names=type, vals=fl )
            rstFile <- file.path(private$projectFolder,paste0(type,".tif"))
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))
        },
        ## split_to_class
        apply_classify = function(layer_name,base_layer,cuts){
            
            ## check base layer exists
            if(!(base_layer %in% names(private$brk))){
                stop(paste(c("Missing layers:",base_layer,sep="\n")))
            }

            ## check layer_name isn't already used
            if(layer_name %in% names(private$brk)){
                stop("layer_name is already used")
            }

            ## Check layer_name isn't reserved
            if(layer_name %in% names(private$reserved_layers)){
                stop("layer_name is reserved")
            }
            
            ## load base layer and mask out channel
            x <-  terra::mask( private$brk[[base_layer]], private$brk[["channel"]], inverse=TRUE)

            ## work out breaks
            brk <- as.numeric(cuts)
            if( length(brk)==1 ){
                ## this defines brks in the same way as cut would otherwise
                rng <- as.numeric( terra::global(x, fun="range",na.rm=TRUE) )
                brk <- seq(rng[1],rng[2],length=brk+1)
            }
            if( any(is.na(brk)) ){ stop("NA value in brk") }

            ## cut the raster
            outFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))

            
            private$brk <- c( private$brk,
                             terra::classify(x,rcl=brk,include.lowest=TRUE,
                                             filename= outFile, names=layer_name))

            out <- list(type="classification", layer=base_layer, cuts=brk)
            writeLines( jsonlite::toJSON(out), file.path(private$projectFolder,paste0(layer_name,".json")) )
        },
        ## split_to_class
        apply_combine_classes = function(layer_name,pairs,burns){
            
            ## check all cuts and burns are in possible layers
            rq <- c(pairs,burns)
            has_rq <- rq %in% names(private$brk)
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }

            ## check layer_name isn't already used
            if(layer_name %in% names(private$brk)){
                stop("layer_name is already used")
            }
                
            ## work out new pairings by cantor method then renumber
            init <- TRUE
            for(ii in pairs){
                x <-  terra::mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
                if(init){
                    cp <- x
                    init <- FALSE
                }else{                    
                    cp <- 0.5*(cp+x)*(cp+x+1)+x
                    uq <- sort(terra::unique(cp)[[1]])
                    cp <- terra::classify(cp, c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf),include.lowest=TRUE)
                }
            }

            ## put all the burns into a single raster
            brn <- cp; brn[] <- NA
            for(ii in burns){
                x <-  terra::mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
                brn <- terra::cover(x,brn)
            }
            ## add burns to pairs
            cp <- terra::cover(brn,cp)
            uq <- sort(terra::unique(cp)[[1]])
            cts <- c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf)
            cp <- terra::classify(cp,cts) + 1 ## classify returns numeric values starting at 0
            if(length(burns)>0){ brn <- terra::classify(brn,cts) +1 }
            
            
            ## TODO - replace with zone taking modal value
            ## make table of layer values - should be able to combine with above??
            cpv <- terra::as.matrix(cp, wide=TRUE) ## quicker when a vector
            uq <- sort(terra::unique(cp)[[1]]) ## unique values
            uqb <- terra::unique(brn)[[1]] ## unique burn values
            
            cuq <- rep(NA,length(uq)) ##index of unique values
            for(ii in which(is.finite(cpv))){
                jj <- cpv[ii]
                if(is.na(cuq[jj])){ cuq[jj] <- ii }
            }
            
            if(!all(is.finite(cuq))){
                stop("Error in computing combinations")
            }

            ## create data frame
            df <- data.frame(uq); names(df) <- layer_name
            for(ii in pairs){
                df[[ii]] <- terra::as.matrix(private$brk[[ii]],wide=TRUE)[cuq] ## read in raster
            }
            df$burns <- df[[layer_name]] %in% uqb

            outFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            private$brk <- c( private$brk, terra::writeRaster( cp, outFile, names=layer_name))
            
            out <- list(type="combination",
                        groups=df)
            writeLines( jsonlite::toJSON(out), file.path(private$projectFolder,paste0(layer_name,".json")) )
            
        },
        
        ## create a model  
        apply_create_model = function(class_lyr,
##                                      dist_lyr,
##                                      delta_dist,
                                      rain_lyr,rainfall_label,
                                      pet_lyr,pet_label,
                                      layer_name,verbose,
                                      sf_opt,
                                      sz_opt){
            ## check layers
            rq <- c("gradient","channel",
                    "band",class_lyr,
                    rain_lyr,pet_lyr)
            has_rq <- rq %in% names(private$brk)            
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }

            if(verbose){ cat("Checking model layer","\n") }
            jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[class_lyr]])),".json")
            if( !file.exists(jsonFile) ){
                stop("No json file giving basis of the classifications")
            }
            json <- jsonlite::fromJSON(jsonFile)
            if( !(json$type=="combination") ){
                stop("Model layer id not a combinations")
            }
            json <- json$group
            if( !("band" %in% names(json)) ){ stop("Model layer should be a combination with the band") }


            if(verbose){ cat("Loading layers","\n") }
            M <- list()
            for(ii in rq){
                M[[ii]] <- terra::as.matrix( private$brk[[ii]], wide=TRUE )
            }
        
            ## check flow routing
            rq <- c( file.path(private$projectFolder,"dem.rds"),
                    file.path(private$projectFolder,"channel.rds") )
            if( ! all( file.exists(rq) ) ){
                stop("No flow routing records defined\n",
                     "Try running compute_flow_paths first")
            }else{
                if(verbose){ cat("Loading routing","\n") }
                flow_routing <- readRDS(rq[1])
                channel_routing <- readRDS( rq[2] )
            }
            
            ## make basic template based on sf_opt and sz_opt
            tmplate <- list(id = integer(0),
                            band = integer(0),
                            states = setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz")),
                            properties = setNames(rep(0,4), c("area","width","Dx","gradient")),
                            sf = list(),
                            rz = list(type="orig", parameters = c("s_rzmax" = 0.1)),
                            uz = list(type="orig", parameters = c("t_d" = 8*60*60)),
                            sz = list(),
                            sf_flow_direction = numeric(0), #list(id = integer(0), fraction = numeric(0)),
                            sz_flow_direction = numeric(0), #list(id = integer(0), fraction = numeric(0)),
                            initialisation = c("s_rz_0" = 0.75, "r_uz_sz_0" = 1e-7),
                            precip = numeric(0),
                            pet = numeric(0),
                            class = list()
                            )
            
            tmplate$sf <- switch(sf_opt,
                                 "cnst" = list(type = "cnst",
                                               parameters = c("c_sf" = 0.3, "d_sf" = 0.0,
                                                              "s_raf" = 0.0, "t_raf" = 999.9)),
                                 "kin" = list(type = "kin",
                                              parameters = c("n" = 0.03,
                                                             "s_raf" = 0.0, "t_raf" = 999.9)),
                                 stop("Unrecognised surface option")
                                 )
            
            tmplate$sz <- switch(sz_opt,
                                 "exp" = list(type = "exp",
                                              parameters = c( "t_0" = 0.135, "m" = 0.04 )),
                                 "bexp" = list(type = "bexp",
                                               parameters = c( "t_0" = 0.135, "m" = 0.04 , "h_sz_max" = 5)),
                                 "dexp" = list(type = "dexp",
                                               parameters = c( "t_0" = 0.135, "m" = 0.04, "m2" = 0.1, "omega"=0.5)),
                                 "cnst" = list(type = "cnst",
                                               parameters = c( "v_sz" = 0.1, "h_sz_max" = 5 )),
                                 stop("Unrecognised saturated zone option")
                                 )
            if(is.null(rain_lyr)){ tmplate$precip <- c("precip"=1) }#list(name="precip", fraction = 1) }
            if(is.null(pet_lyr)){ tmplate$pet <- c("pet"=1) }#list(name = "pet", fraction = 1) }
            
            ## work out some properties of the brick
            cell_area <- prod( terra::res(private$brk) )
            nr <- nrow(M[[1]]); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## number the class layer sequentally from 1
            if(verbose){ cat("Checking class is ordered sequentually","\n") }
            M[[class_lyr]][] <- as.integer( M[[class_lyr]] )
            tmp <- sort(unique(as.vector(M[[class_lyr]])))
            if( all(diff(tmp)==1) ){ ## a seq e.g. from classification
                M[[class_lyr]] <- M[[class_lyr]] + (1-tmp[1])
                json$id <- json[[class_lyr]] + (1-tmp[1])
            }else{
                idx <- match(M[[class_lyr]], tmp)
                M[[class_lyr]][] <- seq_along(tmp)[idx]
                idx <- match( json[[class_lyr]], tmp)
                json$id <- seq_along(tmp)[idx]
            }
            
            ## make the initial hru map
            if(verbose){ cat("make inital HRU map","\n") }
            n_channel <- nrow(private$shp)
            chn_idx <- which( is.finite( M[["channel"]] ) )
            hru_map <- M[[class_lyr]] + n_channel
            hru_map[chn_idx] <- M[["channel"]][chn_idx]
            n_hru <- max(hru_map,na.rm=T)
            json$id <- json$id + n_channel
            
            ## initalise the hrus
            if(verbose){ cat("Initialise the HRUs","\n") }
            hru <- rep(list(tmplate), n_hru)

            ## ############################################
            if(verbose){ cat("Processing channel HRU values","\n") }
            ## Loop channel properties
            class_names <- setdiff(names(private$shp), c("id","band","width","length","slope"))
            for(ii in 1:n_channel){
                ## channel_hru[[ii]]$id <- as.integer( private$shp$id[ii] )
                hru[[ii]]$band <- as.integer( private$shp$band[ii] )
                hru[[ii]]$properties["area"] <- 0
                hru[[ii]]$properties["width"] <- as.numeric( private$shp$width[ii] )
                hru[[ii]]$properties["Dx"] <- as.numeric( private$shp$length[ii] )
                hru[[ii]]$properties["gradient"] <- as.numeric( private$shp$slope[ii] )
                hru[[ii]]$class <- as.list( private$shp[ii,class_names] )
            }
            print(length(hru))
            ## loop the channel routing
            for(ii in 1:nrow(channel_routing)){
                jj <- channel_routing[ii,1]
                kk <- paste(channel_routing[ii,2])
                ff <- channel_routing[ii,3]
                hru[[jj]]$sf_flow_direction[kk] <- ff ## since everything appears once in routing
            }
            print(length(hru))
            ## loop to get outlets
            outlets <- list()
            cnt <- 1
            for(ii in 1:n_channel){
                if( length(hru[[jj]]$sf_flow_direction)==0 ){
                    outlets[[cnt]] <- data.frame(names = paste0("q_sf_",ii),
                                                 id = as.integer( ii ),
                                                 flux = "q_sf", scale = 1.0)
                }
            }
            outlets <- do.call(rbind,outlets)
            print(length(hru))
            ## loop the channel cells
            for(ii in chn_idx){
                jj <- hru_map[ii]
                hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] + cell_area
                if( !is.null(rain_lyr) ){
                    kk <- paste( M[[rain_lyr]][ii] )
                    if(!(kk%in%names( hru[[jj]]$precip))){
                        hru[[jj]]$precip[kk] <- 0
                    }
                    hru[[jj]]$precip[kk] <- hru[[jj]]$precip[kk] + 1
                }
                if( !is.null(pet_lyr) ){
                    kk <- paste( M[[pet_lyr]][ii] )
                    if(!(kk%in%names( hru[[jj]]$pet))){
                        hru[[jj]]$pet[kk] <- 0
                    }
                    hru[[jj]]$pet[kk] <- hru[[jj]]$pet[kk] + 1
                }
            }

            ## ############################################
            if(verbose){ cat("Processing hillslope HRU values","\n") }
            print(length(hru))
            ## pass one to do some properties and rainfall
            for(rw in 1:nrow(flow_routing)){
                ii <- flow_routing[rw,1] ## cell number
                jj <- hru_map[ii] ## current hru

                hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] + 1
                hru[[jj]]$properties["gradient"] <- M[["gradient"]][ii]
                hru[[jj]]$band <- min( hru[[jj]]$band , M[["band"]][ii] )

                if( !is.null(rain_lyr) ){
                    kk <- paste( M[[rain_lyr]][ii] )
                    if(!(kk%in%names( hru[[jj]]$precip))){
                        hru[[jj]]$precip[kk] <- 0
                    }
                    hru[[jj]]$precip[kk] <- hru[[jj]]$precip[kk] + 1
                }
                if( !is.null(pet_lyr) ){
                    kk <- paste( M[[pet_lyr]][ii] )
                    if(!(kk%in%names( hru[[jj]]$pet))){
                        hru[[jj]]$pet[kk] <- 0
                    }
                    hru[[jj]]$pet[kk] <- hru[[jj]]$pet[kk] + 1
                }
            }
            ## second pass to do flow routing
            for(rw in 1:nrow(flow_routing)){
                ii <- flow_routing[rw,1] ## cell number
                jj <- hru_map[ii] ## current hru

                if( M[["band"]][ii] > hru[[jj]]$band ){ next } ## look only at lowest band

                gcl <- flow_routing[rw,2]
                dcl <- flow_routing[rw,3]
                w <- flow_routing[rw,4:11]
                
                if( !any(w>0) ){ stop("Hillslope cell with no outflow!") }

                hru[[jj]]$properties["width"] <- hru[[jj]]$properties["width"] + dcl
               
                is_ds <- w>0
                kk <- (ii + delta)[ is_ds ]
                ngh <- paste( hru_map[ kk ] )
                new_ngh <- setdiff( ngh, names( hru[[jj]]$sf_flow_direction ) )
                hru[[jj]]$sf_flow_direction[ new_ngh ] <- 0
                hru[[jj]]$sf_flow_direction[ ngh ] <- hru[[jj]]$sf_flow_direction[ ngh ] + gcl*w[is_ds]
            }
            ## pass over hillslope HRU to tidy up properties
            for(jj in (n_channel+1):n_hru){
                hru[[jj]]$properties["gradient"] <- hru[[jj]]$properties["gradient"] / hru[[jj]]$properties["area"]
                hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] * cell_area
                hru[[jj]]$properties["length"] <- hru[[jj]]$properties["area"] / hru[[jj]]$properties["width"]
                if( length(hru[[jj]]$sf_flow_direction)==0 ){
                    stop("Hillslope HRU with no outflow")
                }
                hru[[jj]]$sf_flow_direction <- hru[[jj]]$sf_flow_direction / sum(hru[[jj]]$sf_flow_direction)
                hru[[jj]]$class <- as.list( json[json$id==jj,] )
                hru[[jj]]$class$id <- NULL
            }
                
            ## ############################################
            if(verbose){ cat("Finalising HRUs","\n") }
            bnd <- rep(NA,n_hru)
            for(ii in 1:n_hru){
                bnd[ii] <- hru[[ii]]$band
                hru[[ii]]$id <- ii
                hru[[ii]]$precip <- list(name = paste0(rainfall_label,names(hru[[ii]]$precip)),
                                         fraction = as.numeric( hru[[ii]]$precip / sum(hru[[ii]]$precip) ))
                hru[[ii]]$pet <- list(name = paste0(pet_label,names(hru[[ii]]$pet)),
                                      fraction = as.numeric( hru[[ii]]$pet / sum(hru[[ii]]$pet) ))
                hru[[ii]]$sf_flow_direction <- list(
                    id = as.integer( names(hru[[ii]]$sf_flow_direction) ),
                    fraction = as.numeric( hru[[ii]]$sf_flow_direction ))
                if( ii > n_channel ){
                    hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
                }else{
                    hru[[ii]]$sz_flow_direction <- list(id = integer(0),fraction=numeric(0))
                }
            }

            ## ###########################################
            ## order the hru's
            ## for historic reasons the HRU id's must be in computational order...
            idx <- order(bnd)
            browser()
            hru_map <- idx[ hru_map ]
            hru <- hru[idx]
            for(ii in 1:n_hru){
               hru[[ii]]$id <- ii
               hru[[ii]]$sf_flow_direction$id <- as.integer( idx[hru[[ii]]$sf_flow_direction$id] )
               hru[[ii]]$sz_flow_direction$id <- as.integer( idx[hru[[ii]]$sz_flow_direction$id] )
            }

            ## make output
            if(verbose){ cat("Making output","\n") }
            rst <- terra::rast( private$brk[[ "channel" ]],
                               names = layer_name, vals=hru_map)
            model <- list(hru=hru, output_flux = outlets)
            model$map <- paste0(layer_name,".tif")
            terra::writeRaster(rst,model$map,overwrite=TRUE)
            saveRDS(model,paste0(layer_name,".rds"))
        }
    )
)


            

                
##                 idx <- channel_routing[,"from"]== private$shp$id[ii]
##                 if( any(idx) ){
##                     channel_hru[[ii]]$sf_flow_direction$id <- as.integer( channel_routing[idx,"to"] )
##                     channel_hru[[ii]]$sf_flow_direction$fraction <- as.integer( channel_routing[idx,"frac"] )
##                     channel_hru[[ii]]$sf_flow_direction <- channel_hru[[ii]]$sz_flow_direction
##                 }
##                 channel_hru[[ii]]$class <- as.list( private$shp[ii,class_names] )

##                 cnt <- cnt + 1
##             }
            
            
                                                                         


            
##             ## process hilldlope HRU's
##             ## read in classification and distance layers
##             cls <-  private$brk[[class_lyr]] ## channel values are NA
##             dst <- private$brk[[dist_lyr]] ## !!!!We assume that all channel pixels have a distance that puts them in the correct order!!!!!
            
##             cls <- cls + max(private$shp$id) ## alter class so greater then river channel id
##             hmap <-  terra::cover( private$brk[["channel"]], cls) ## make map of HRUs - but numbering not yet correct
##             names(hmap) <- "hru"

##             ## work out the order of the hrus and add class information
##             tbl <- terra::zonal(dst,hmap,min) # minimum distance for each hillslope classification
##             tbl$min_dst <- tbl[[dist_lyr]]
##             tbl[[dist_lyr]] <- NULL
   
##             ## add any missing channels
##             if( !all(private$shp$id %in% tbl$hru ) ){ stop("where is the HRU!!!") }

##             ## add class information
##             jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[class_lyr]])),".json")
##             if( !file.exists(jsonFile) ){
##                 warning("No json file giving basis of the classifications")
##             }else{
##                 tmp <- jsonlite::fromJSON(jsonFile)$groups
##                 tmp$hru <- tmp[[class_lyr]] + max(private$shp$id)
##                 tbl <- merge(tbl, tmp, by="hru", all.x=TRUE)
##             }
            
##             ##tbl$is_channel <- tbl$hru %in% private$shp$id
##             tbl <- merge(tbl,as.data.frame(private$shp),by.x="hru",by.y="id",all.x=TRUE) ## add channel information

##             ## order and renumber
##             tbl <- tbl[ order(tbl$min_dst), ]
##             hmap <- terra::subst(hmap, tbl$hru, 0:(nrow(tbl)-1)) ## this line is slow
##             tbl$hru <- 0:(nrow(tbl)-1)

##             ## make hrus
##             hru <- lapply(1:nrow(tbl), function(ii){ ## slow but not as slow as subst..
##                 tmp <- as.list(tbl[ii,])
##                 out <- tmplate
##                 out$id <- tmp$hru
##                 tmp$hru <- NULL
##                 out$class <- tmp
##                 return(out)
##             })
            
## ##            })) ## end of first system.time

## ##            print(system.time({
                
##             ## work out the properties
##             if(verbose){ cat("Computing properties","\n") }
     
##             ## it is quickest to compute using blocks of a raster
##             ## however for small rasters we will just treat as a single block
##             d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE )
##             mp <- terra::as.matrix( hmap, wide=TRUE )
##             dst <- terra::as.matrix( dst, wide=TRUE )
##             gr <- terra::as.matrix( private$brk[["gradient"]], wide=TRUE )
##             ctch <- terra::as.matrix( private$brk[["catchment"]],  wide=TRUE )
##             if( !is.null(rain_lyr) ){
##                 rain <- terra::as.matrix( private$brk[[rain_lyr]], wide=TRUE )
##             }
##             if( !is.null(pet_lyr) ){
##                 pet <- terra::as.matrix( private$brk[[pet_lyr]], wide=TRUE )
##             }

## ##            })) ## end of second system.time
##             ##atb <- terra::as.matrix( private$brk[["atb"]]  , wide=TRUE )

## ##            print(system.time({




            
##             if( verbose ){ print("Computing hillslope routing") }
##             ## distances and contour lengths
##             ## distance between cell centres
##             rs <- terra::res( private$brk )
##             dxy <- rep(sqrt(sum(rs^2)),8)
##             dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
##             dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
##             nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
##             ## if we go up in height order then we are working from near the channel to the heighest point
##             idx <- order(d,na.last=NA)
##             n_to_eval <- length(idx)
##             cellArea <- prod(rs)
            
##             for(ii in idx){
##                 id <- mp[ii]
##                 jj <- id+1 ## index in the hru list (id +1)
                
##                 hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] + 1
##                 ##hru[[jj]]$properties["atb_bar"] <- hru[[jj]]$properties["atb_bar"] + atb[ii]
##                 hru[[jj]]$properties["gradient"] <- hru[[jj]]$properties["gradient"] + gr[ii]
                
##                 ## work out subsurface flow for non-channel
##                 if( is.na(hru[[jj]]$class$startNode) & (dst[ii] <= (hru[[jj]]$class$min_dst + delta_dist)) ){
                    
##                     ## look for flow direction
##                     ngh <- ii + delta ## neighbouring cells
##                     nghid <- mp[ngh]
##                     cjdx <- ctch[ngh]
                    
##                     ## compute gradient
##                     grd <- (d[ii]-d[ngh])/dxy
##                     ## test if can be used
##                     to_use <- is.finite(grd) & (grd > 0) & is.finite(nghid) & (nghid < id) & (cjdx==ctch[ii])
                    
##                     if( any(to_use) ){
##                         nghid <- paste(nghid[to_use])
##                         tmp <- nghid[!(nghid %in% names(hru[[jj]]$sz_flow_direction))] 
##                         if(length(tmp)>0){
##                             hru[[jj]]$sz_flow_direction[ tmp ] <- 0
##                         }
##                         tmp <- tapply(grd[to_use]*dcl[to_use],nghid,sum)
##                         hru[[jj]]$sz_flow_direction[ names(tmp) ] <-
##                             hru[[jj]]$sz_flow_direction[ names(tmp) ] + tmp
##                         hru[[jj]]$properties["width"] <- hru[[jj]]$properties["width"] + sum(dcl[to_use])
##                     }
                    
##                 }

##                 ## work out precipitation
##                 if(!is.null(rain_lyr) && is.finite(rain[ii])){
##                     nm <- paste0(rainfall_label,rain[ii])
##                     if(!(nm %in% names(hru[[jj]]$precip))){
##                         hru[[jj]]$precip[nm] <- 0
##                     }
##                     hru[[jj]]$precip[nm] <- hru[[jj]]$precip[nm] + 1
##                 }
                
##                 ## work out pet
##                 if(!is.null(pet_lyr) && is.finite(pet[ii])){
##                     nm <- paste0(pet_label,pet[ii])
##                     if(!(nm %in% names(hru[[jj]]$pet))){ 
##                         hru[[jj]]$pet[nm] <- 0
##                     }
##                     hru[[jj]]$pet[nm] <- hru[[jj]]$pet[nm] + 1
##                 }   
                    
##                 n_to_eval[2] <- n_to_eval[2] + 1
##                 if( verbose & (n_to_eval[2] > n_to_eval[3]) ){
##                     cat(round(100*n_to_eval[3] / n_to_eval[1],1),
##                         "% complete","\n")
##                     n_to_eval[3] <- n_to_eval[3] + n_to_eval[1]/20
##                 }
                
##             }

## ##            })) ## end of third system.time

## ##            print(system.time({
                
##             ## second pass to correct sumations and compute surface
##             sN <- sapply(hru,function(h){h$class$startNode})

##             outlet_id <- NULL
##             no_outflow <- NULL
##             for(ii in 1:length(hru)){

##                 ## check precip
##                 if(!is.null(rain_lyr)){
##                     if( sum(hru[[ii]]$precip) != hru[[ii]]$properties["area"] ){
##                         warning(paste("HRU",hru[[ii]]$id," - Precip cell count is",sum(hru[[ii]]$precip),
##                                       "full cell count is",hru[[ii]]$properties["area"]))
##                     }
##                 }
##                 hru[[ii]]$precip <- list(name = as.character(names(hru[[ii]]$precip)),
##                                          fraction = as.numeric(hru[[ii]]$precip/sum(hru[[ii]]$precip)) )
                
##                 ## check pet
##                 if(!is.null(pet_lyr)){
##                     if( sum(hru[[ii]]$pet) != hru[[ii]]$properties["area"] ){
##                         warning(paste("HRU",hru[[ii]]$id," - PET cell count is",sum(hru[[ii]]$precip),
##                                       "full cell count is",hru[[ii]]$properties["area"]))
##                     }
##                 }
                
##                 hru[[ii]]$pet <- list(name = as.character(names(hru[[ii]]$pet)),
##                                       fraction = as.numeric(hru[[ii]]$pet/sum(hru[[ii]]$pet)) )
                
##                 ##hru[[ii]]$properties["atb_bar"] <- hru[[ii]]$properties["atb_bar"] / hru[[ii]]$properties["area"]
##                 hru[[ii]]$properties["gradient"] <- hru[[ii]]$properties["gradient"] / hru[[ii]]$properties["area"]
##                 hru[[ii]]$properties["area"] <- hru[[ii]]$properties["area"] * cellArea
                
##                 if( !is.na(hru[[ii]]$class$startNode) ){
##                     ## then a channel
##                     hru[[ii]]$properties["Dx"] <- hru[[ii]]$class$length
##                     hru[[ii]]$properties["width"] <- hru[[ii]]$class$width
##                     hru[[ii]]$properties["gradient"] <- hru[[ii]]$class$slope ## PJS TODO document this and requirement to have it as a variable
##                     jdx <- which( sN == hru[[ii]]$class$endNode)
##                     if(length(jdx) >0){
##                         ## has downstream connenction
##                         hru[[ii]]$sf_flow_direction <- list(
##                             id = as.integer( jdx - 1 ),
##                             fraction = rep( as.numeric( 1 / length(jdx) ),length(jdx) ))
##                     }else{
##                         ## goes to an outlet
##                         hru[[ii]]$sf_flow_direction <- list(id = integer(0),fraction=numeric(0))
##                         outlet_id <- c(outlet_id,hru[[ii]]$id)
##                     }
##                     ## set subsruface to match surface
##                     hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
                    
##                 }else{
##                     if( length(hru[[ii]]$sz_flow_direction)==0 ){
##                         no_outflow <- c(no_outflow, hru[[ii]]$id)
##                     }
                    
##                     hru[[ii]]$sz_flow_direction <- list(
##                         id = as.integer(names( hru[[ii]]$sz_flow_direction )),
##                         fraction = as.numeric( hru[[ii]]$sz_flow_direction/sum(hru[[ii]]$sz_flow_direction) ))

##                     ## set surface to match subsurface
##                     hru[[ii]]$sf_flow_direction <- hru[[ii]]$sz_flow_direction
##                     hru[[ii]]$properties["Dx"] <- hru[[ii]]$properties["area"]/hru[[ii]]$properties["width"]
##                 }
##             }

## ##            })) ## end of system.time
            
##             ## stop if any points with no outflow that aren't channels
##             if( length(no_outflow) >0 ){
##                 stop( "The following HRUs have no valid outflows and are not channels: \n",
##                      paste(no_outflow, collapse=", "))
##             }
            
##             ## ############################################
##             ## Add outlet at all outlets from river network
##             ## ############################################
##             ##idx <- sapply(hru,function(x){ifelse(length(x$sf_flow_direction$id)==0,x$id,NULL)})
##             output_flux<- data.frame(name = paste0("q_sf_",outlet_id),
##                                            id = as.integer( outlet_id ), flux = "q_sf", scale = 1.0)
                        
##             ## ############################
##             ## save model
##             model <- list(hru=hru, output_flux = output_flux)
##             model$map <- paste0(layer_name,".tif")
##             terra::writeRaster(hmap,model$map,overwrite=TRUE)
##             saveRDS(model,paste0(layer_name,".rds"))
##         },
##         ## create a gauge list
##         apply_locate_gauges = function(gauge_map,
##                                        gauge_name
##                                        ){
##             if( length(private$shp) == 0 ){
##                 stop("Channel needs to be added")
##             }
            
##             dst <- distance(gauge_map,private$shp) ## distance between gauges(row) and channels(column)
            
##             idx <- apply(dst,1,which.min) ## index of closest channel for each gauge
##             gauge_vec$chn_identifier <- chn_poly$name[ apply(dst,1,which.min) ]
##             gauge_vec$chn_distance <- apply(dst,1,min) ## distance to the channel

##         }
##     )
##     )

##     ##     ## #######################
##     ##     ## create a model version 2
##     ##     apply_create_model_new = function(class_lyr,
##     ##                                       rain_lyr,rainfall_label,
##     ##                                       pet_lyr,pet_label,
##     ##                                       layer_name,verbose,
##     ##                                       sf_opt,
##     ##                                       sz_opt){
    
    
##     ##         rq <- c("gradient","channel","band",
##     ##                 class_lyr,
##     ##                 rain_lyr,pet_lyr)
##     ##         has_rq <- rq %in% names(private$brk)            
##     ##         if(!all(has_rq)){
##     ##             stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
##     ##         }
    
##     ##         rq <- c( file.path(private$projectFolder,"flow_direction.rds"),
##     ##                 file.path(private$projectFolder,"channel_direction.rds") )
##     ##         if( ! all( file.exists(rq) ) ){
##     ##             stop("No flow routing records defined\n",
##     ##                  "Try running compute_flow_paths first")
##     ##         }
    
##     ##         if(verbose){ cat("Creating HRU map","\n") }
##     ##         nchn <- max(private$shp$id)
##     ##         ## make map of hrus
##     ##         hmap <-  private$brk[[class_lyr]]
##     ##         hmap <- hmap + nchn ## alter class so greater then river channel id
##     ##         hmap <-  terra::cover( private$brk[["channel"]], hmap) ## make map of HRUs - but numbering not yet correct
##     ##         names(hmap) <- "hru"
##     ##         nhru <- terra::global(hmap,max,na.rm=T)
    
##     ##         if(verbose){ cat("Setting up HRUs","\n") }
        
##     ##         ## make basic template based on sf_opt and sz_opt
##     ##         tmp_sf <- switch(sf_opt,
##     ##                          "cnst" = list(type = "cnst",
##     ##                                        parameters = c("c_sf" = 0.3, "d_sf" = 0.0,
##     ##                                                       "s_raf" = 0.0, "t_raf" = 999.9)),
##     ##                          "kin" = list(type = "kin",
##     ##                                       parameters = c("n" = 0.03,
##     ##                                                      "s_raf" = 0.0, "t_raf" = 999.9)),
##     ##                          stop("Unrecognised surface option")
##     ##                          )
##     ##         tmp_sz <- switch(sz_opt,
##     ##                          "exp" = list(type = "exp",
##     ##                                       parameters = c( "t_0" = 0.135, "m" = 0.04 )),
##     ##                          "bexp" = list(type = "bexp",
##     ##                                        parameters = c( "t_0" = 0.135, "m" = 0.04 , "h_sz_max" = 5)),
##     ##                          "dexp" = list(type = "dexp",
##     ##                                        parameters = c( "t_0" = 0.135, "m" = 0.04, "m2" = 0.1, "omega"=0.5)),
##     ##                          "cnst" = list(type = "cnst",
##     ##                                        parameters = c( "v_sz" = 0.1, "h_sz_max" = 5 )),
##     ##                          stop("Unrecognised saturated zone option")
##     ##                          )
##     ##         if(is.null(rain_lyr)){ tmp_prcp <- c("precip"=0) }else{ tmp_prcp <- numeric(0) }
##     ##         if(is.null(pet_lyr)){ tmp_pet <- c("pet"=0) }else{ tmp_pet <- numeric(0) }
##     ##         tmplate <- list(id = as.integer(-999),
##     ##                         band = as.integer(Inf),
##     ##                         states = setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz")),
##     ##                         properties = setNames(rep(0,4), c("area","width","Dx","gradient")),
##     ##                         sf = tmp_sf,
##     ##                         rz = list(type="orig", parameters = c("s_rzmax" = 0.1)),
##     ##                         uz = list(type="orig", parameters = c("t_d" = 8*60*60)),
##     ##                         sz = tmp_sz,
##     ##                         sf_flow_direction = numeric(0), #list(id = integer(0), fraction = numeric(0)),
##     ##                         sz_flow_direction = numeric(0), #list(id = integer(0), fraction = numeric(0)),
##     ##                         initialisation = c("s_rz_0" = 0.75, "r_uz_sz_0" = 1e-7),
##     ##                         precip = tmp_prcp,
##     ##                         pet = tmp_pet,
##     ##                         class = list()
##     ##                         )
##     ##         hru <- rep(list(tmplate),nhru)
    
##     ##         ## open required maps
##     ##         mp <- terra::as.matrix( hmap, wide=TRUE )
##     ##         gr <- terra::as.matrix( private$brk[["gradient"]], wide=TRUE )
##     ##         bnd <- terra::as.matrix( private$brk[["band"]], wide=TRUE )
##     ##         if( !is.null(rain_lyr) ){
##     ##             rain <- terra::as.matrix( private$brk[[rain_lyr]], wide=TRUE )
##     ##         }
##     ##         if( !is.null(pet_lyr) ){
##     ##             pet <- terra::as.matrix( private$brk[[pet_lyr]], wide=TRUE )
##     ##         }
    
##     ##         ## compute areas, gradients and rainfalls
##     ##         if(verbose){ print("Computing summaries") }
    
##     ##         idx <-  which(is.finite(mp))
##     ##         n_to_eval <- c(length(idx),0,length(idx)/20)
##     ##         for(ii in idx){
##     ##             id <- mp[ii]
    
##     ##             hru[[id]]$properties["area"] <- hru[[id]]$properties["area"] + 1
##     ##             hru[[id]]$properties["gradient"] <- hru[[id]]$properties["gradient"] + grad[ii]
##     ##             if( hru[[id]]$band > bnd[ii] ){ hru[[id]]$band <- bnd[ii] }
    
    
##     ##             ## work out precipitation
##     ##             if(!is.null(rain_lyr)){ 
##     ##                 nm <- paste0(rainfall_label,rain[ii])
##     ##                 if(!(nm %in% names(hru[[jj]]$precip))){ hru[[id]]$precip[nm] <- 0 }
##     ##             }else{
##     ##                 nm <- "precip"
##     ##             }
##     ##             hru[[id]]$precip[nm] <- hru[[id]]$precip[nm] + 1
    
    
##     ##             ## work out pet
##     ##             if(!is.null(pet_lyr)){
##     ##                 nm <- paste0(pet_label,pet[ii])
##     ##                 if(!(nm %in% names(hru[[jj]]$pet))){ hru[[id]]$pet[nm] <- 0 }
##     ##             }else{
##     ##                 nm <- "pet"
##     ##             }
##     ##             hru[[id]]$pet[nm] <- hru[[id]]$pet[nm] + 1
    
    
##     ##             n_to_eval[2] <- n_to_eval[2] + 1
##     ##             if( verbose & (n_to_eval[2] > n_to_eval[3]) ){
##     ##                 cat(round(100*n_to_eval[3] / n_to_eval[1],1),
##     ##                     "% complete","\n")
##     ##                 n_to_eval[3] <- n_to_eval[3] + n_to_eval[1]/20
##     ##             }
    
##     ##         }
    
##     ##         ## do hillslope flow routing
##     ##         if(verbose){ cat("Computing hillslope flow directions and inputs","\n") }
    
##     ##         ## flow direction
##     ##         fd <- file.path(private$projectFolder,"flow_direction.rds")
    
##     ##         ## distances and contour lengths
##     ##         rs <- terra::res( private$brk )
##     ##                                         #dxy <- rep(sqrt(sum(rs^2)),8)
##     ##                                         #dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
##     ##         dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
##     ##         nr <- nrow(mp); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
##     ##         cellArea <- prod(rs)
    
##     ##         n_to_eval <- c(nrow(fd),0,nrow(fd)/20)
    
##     ##         for(rw in nrow(fd):1){
##     ##             ii <- fd[rw,1]
##     ##             w <- fd[rw,-1]
##     ##             id <- mp[ii]
    
##     ##             ## flow direction
##     ##             jj <- ii+delta
##     ##             to_use <- (bnd[jj] < bnd[ii]) & (w>0)
    
##     ##             if( bnd[ii] == hru[[id]]$band ){
    
##     ##                 jj <- ii+delta
##     ##                 to_use <- (bnd[jj] < bnd[ii]) & (w>0)
##     ##                 if( any(to_use) ){
    
##     ##                     jj <- jj[to_use]
##     ##                     cl <- dcl[to_use]
##     ##                     wcl <- cl*w[to_use]
##     ##                     did <- mp[jj]
##     ##                     swcl <- tapply(wcl,did,sum)
    
##     ##                     idx <- setdiff(names(swcl),names( hru[[id]]$sf_flow_direction ))
##     ##                     hru[[id]]$sf_flow_direction[ idx ] <- 0
##     ##                     hru[[id]]$sf_flow_direction[ names(swcl) ] <- hru[[id]]$sf_flow_direction[ names(swcl) ] + swcl
##     ##                     hru[[id]]$properties["width"] <- hru[[id]]$properties["width"] + sum(cl)
##     ##                 }
##     ##             }
    
##     ##             n_to_eval[2] <- n_to_eval[2] + 1
##     ##             if( verbose & (n_to_eval[2] > n_to_eval[3]) ){
##     ##                 cat(round(100*n_to_eval[3] / n_to_eval[1],1),
##     ##                     "% complete","\n")
##     ##                 n_to_eval[3] <- n_to_eval[3] + n_to_eval[1]/20
##     ##             }
    
##     ##         }
    
##     ##         ## Add class infomation
##     ##         if(verbose){ print("Adding class infomation and standardising") }
##     ##         cls_data <- self$get_method(cls_lyr)$groups
##     ##         chn_data <- as.data.frame( private$shp )
##     ##         for(ii in 1:nhru){
##     ##             hru[[ii]]$id <- ii
    
##     ##             n <- sum(hru[[ii]]$precip)
##     ##             if(n!=hru[[ii]]$properties["area"]){
##     ##                 stop("Incorrect rainfall area")
##     ##             }
##     ##             hru[[ii]]$precip <- hru[[ii]]$precip / hru[[ii]]$properties["area"]
##     ##             hru[[ii]]$precip <- list(name = as.character(names(hru[[ii]]$precip)),
##     ##                                      fraction = as.numeric(hru[[ii]]$precip) )
    
##     ##             n <- sum(hru[[ii]]$pet)
##     ##             if(n!=hru[[ii]]$properties["area"]){
##     ##                 stop("Incorrect PETl area")
##     ##             }
##     ##             hru[[ii]]$pet <- hru[[ii]]$pet / hru[[ii]]$properties["area"]
##     ##             hru[[ii]]$pet <- list(name = as.character(names(hru[[ii]]$pet)),
##     ##                                   fraction = as.numeric(hru[[ii]]$pet) )
    
##     ##             hru[[ii]]$properties["gradient"] <- hru[[ii]]$properties["gradient"] / hru[[ii]]$properties["area"]
    
##     ##             hru[[ii]]$properties["area"] <- hru[[ii]]$properties["area"]*cellArea
    
##     ##             if(ii > nchn){
##     ##                 hru[[ii]]$class <- as.list( cls_data[cls_data[[cls_lyr]]==ii-nchn,] )
##     ##                 hru[[ii]]$sf_flow_direction <- hru[[ii]]$sf_flow_direction / sum( hru[[ii]]$sf_flow_direction )
##     ##                 hru[[ii]]$sf_flow_direction <- list(id = as.integer(
##     ##                 hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
##     ##             }else{
##     ##                 hru[[ii]]$class <- as.list( chn_data[chn_data$id ==ii,] )
##     ##                 hru[[ii]]$properties["gradient"] <- hru[[ii]]$class$slope
##     ##                 hru[[ii]]$properties["width"] <- hru[[ii]]$class$width
##     ##                 hru[[ii]]$properties["Dx"] <- hru[[ii]]$class$length
##     ##             }
##     ##         }
    
##     ##         ## pass through second pass to correct sumations and compute surface
##     ##         sN <- sapply(hru,function(h){h$class$startNode})
    
##     ##             outlet_id <- NULL
##     ##             no_outflow <- NULL
##     ##             for(ii in 1:length(hru)){

##     ##                 ## check precip
##     ##                 if(!is.null(rain_lyr)){
##     ##                     if( sum(hru[[ii]]$precip) != hru[[ii]]$properties["area"] ){
##     ##                         warning(paste("HRU",hru[[ii]]$id," - Precip cell count is",sum(hru[[ii]]$precip),
##     ##                                       "full cell count is",hru[[ii]]$properties["area"]))
##     ##                     }
##     ##                 }
##     ##                 hru[[ii]]$precip <- list(name = as.character(names(hru[[ii]]$precip)),
##     ##                                          fraction = as.numeric(hru[[ii]]$precip/sum(hru[[ii]]$precip)) )
    
##     ##                 ## check pet
##     ##                 if(!is.null(pet_lyr)){
##     ##                     if( sum(hru[[ii]]$pet) != hru[[ii]]$properties["area"] ){
##     ##                         warning(paste("HRU",hru[[ii]]$id," - PET cell count is",sum(hru[[ii]]$precip),
##     ##                                       "full cell count is",hru[[ii]]$properties["area"]))
##     ##                     }
##     ##                 }
    
##     ##                 hru[[ii]]$pet <- list(name = as.character(names(hru[[ii]]$pet)),
##     ##                                       fraction = as.numeric(hru[[ii]]$pet/sum(hru[[ii]]$pet)) )
                
## ##                 ##hru[[ii]]$properties["atb_bar"] <- hru[[ii]]$properties["atb_bar"] / hru[[ii]]$properties["area"]
## ##                 hru[[ii]]$properties["gradient"] <- hru[[ii]]$properties["gradient"] / hru[[ii]]$properties["area"]
## ##                 hru[[ii]]$properties["area"] <- hru[[ii]]$properties["area"] * cellArea
                
## ##                 if( !is.na(hru[[ii]]$class$startNode) ){
## ##                     ## then a channel
## ##                     hru[[ii]]$properties["Dx"] <- hru[[ii]]$class$length
## ##                     hru[[ii]]$properties["width"] <- hru[[ii]]$class$width
## ##                     hru[[ii]]$properties["gradient"] <- hru[[ii]]$class$slope ## PJS TODO document this and requirement to have it as a variable
## ##                     jdx <- which( sN == hru[[ii]]$class$endNode)
## ##                     if(length(jdx) >0){
## ##                         ## has downstream connenction
## ##                         hru[[ii]]$sf_flow_direction <- list(
## ##                             id = as.integer( jdx - 1 ),
## ##                             fraction = rep( as.numeric( 1 / length(jdx) ),length(jdx) ))
## ##                     }else{
## ##                         ## goes to an outlet
## ##                         hru[[ii]]$sf_flow_direction <- list(id = integer(0),fraction=numeric(0))
## ##                         outlet_id <- c(outlet_id,hru[[ii]]$id)
## ##                     }
## ##                     ## set subsruface to match surface
## ##                     hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
                    
## ##                 }else{
## ##                     if( length(hru[[ii]]$sz_flow_direction)==0 ){
## ##                         no_outflow <- c(no_outflow, hru[[ii]]$id)
## ##                     }
                    
## ##                     hru[[ii]]$sz_flow_direction <- list(
## ##                         id = as.integer(names( hru[[ii]]$sz_flow_direction )),
## ##                         fraction = as.numeric( hru[[ii]]$sz_flow_direction/sum(hru[[ii]]$sz_flow_direction) ))

## ##                     ## set surface to match subsurface
## ##                     hru[[ii]]$sf_flow_direction <- hru[[ii]]$sz_flow_direction
## ##                     hru[[ii]]$properties["Dx"] <- hru[[ii]]$properties["area"]/hru[[ii]]$properties["width"]
## ##                 }
## ##             }

## ## ##            })) ## end of system.time
            
## ##             ## stop if any points with no outflow that aren't channels
## ##             if( length(no_outflow) >0 ){
## ##                 stop( "The following HRUs have no valid outflows and are not channels: \n",
## ##                      paste(no_outflow, collapse=", "))
## ##             }
            
## ##             ## ############################################
## ##             ## Add outlet at all outlets from river network
## ##             ## ############################################
## ##             ##idx <- sapply(hru,function(x){ifelse(length(x$sf_flow_direction$id)==0,x$id,NULL)})
## ##             output_flux<- data.frame(name = paste0("q_sf_",outlet_id),
## ##                                            id = as.integer( outlet_id ), flux = "q_sf")
                        
## ##             ## ############################
## ##             ## save model
## ##             model <- list(hru=hru, output_flux = output_flux)
## ##             model$map <- paste0(layer_name,".tif")
## ##             terra::writeRaster(hmap,model$map,overwrite=TRUE)
## ##             saveRDS(model,paste0(layer_name,".rds"))
##     ##         }
    
## ##    )

    
