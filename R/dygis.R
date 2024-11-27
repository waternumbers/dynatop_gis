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
#' ctch$compute_band()
#' ctch$compute_properties() # like topograpihc index and contour length
#' ctch$compute_flow_lengths()
#' }
#' ## classify and create a model
#' \donttest{
#' ctch$classify("atb_20","atb",cuts=20) # classify using the topographic index
#' ctch$get_method("atb_20") ## see the details of the classification
#' ctch$combine_classes("atb_20_band",c("atb_20","band")) ## combine classes
#' ctch$create_model(file.path(demo_dir,"new_model"),"atb_20") ## create a model
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
        version = "0.3.2",
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
            ## mask layer to catchment
            layer <- terra::mask(layer,private$brk[[ "catchment" ]])
            ## check no NA values
            ##stopifnot( "Layer should have no NA values within the catchment" =
            ##               terra::global(is.na(layer) & !is.na(private$brk[["catchment"]]),max)==0 )
            
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
                grd <- (d[ii]-d[jdx])/dxy
                gcl <- grd*dcl
                cjdx <- ctch[jdx]
                is_lower <- is.finite(gcl) & gcl>0 & is.finite(cjdx) & cjdx==ctch[ii]
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

            
            ## cut the raster and save
            outFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            private$brk <- c( private$brk,
                             terra::classify(x,rcl=brk,include.lowest=TRUE,
                                             filename= outFile, names=layer_name))
            ## create the output json table
            uq <- unlist(terra::global( private$brk[[layer_name]], range, na.rm=T ))
            ##uq <- sort(terra::unique( private$brk[[layer_name]] )) ## unique values
            out <- list( groups = data.frame( uq[1]:uq[2] ),
                        cuts = list( list(layer = base_layer,cuts=brk) ))
            names(out$groups) <- layer_name
            names(out$cuts) <- layer_name
            ## out <- list(type="classification", layer=base_layer, cuts=brk)
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
            
            out <- list(groups=df)
            writeLines( jsonlite::toJSON(out), file.path(private$projectFolder,paste0(layer_name,".json")) )
            
        },
        
        ## create a model  
        apply_create_model = function(class_lyr,
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
            json <- json$groups

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
            
            ## work out some properties of the brick and channel
            cell_area <- prod( terra::res(private$brk) )
            nr <- nrow(M[[1]]); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            n_channel <- nrow(private$shp)
            
            ## make HRU id map by crossing the class layer and band
            hru_map <- M[[class_lyr]]
            hru_map <- 0.5*(hru_map + M[["band"]])*(hru_map + M[["band"]]+1) + M[["band"]]
            ##hru_map <- M[[class_lyr]]
            tmp <- tapply(M[["band"]],hru_map,max)
            tmp <- sort(tmp)
            tmp <- as.integer(names(tmp)) ## this is original cantor value numbers sorted by band
            id <- ((n_channel+1):(n_channel+length(tmp))) ## corresponding new id
            idx <- match(hru_map,tmp)
            hru_map[] <- id[idx]
            hru_map <- pmin( M[["channel"]], hru_map, na.rm=TRUE )
            
            n_hru <- max(id) ## larget hru id number
            
            ## initalise the hrus
            if(verbose){ cat("Initialise the HRUs","\n") }
            hru <- rep(list(tmplate), n_hru)

            ## ############################################
            if(verbose){ cat("Processing channel HRU values","\n") }
            ## Loop channel properties
            class_names <- setdiff(names(private$shp), c("id","band","width","length","slope"))
            for(ii in 1:n_channel){
                hru[[ii]]$id <- as.integer( private$shp$id[ii] )
                hru[[ii]]$band <- as.integer( private$shp$band[ii] )
                hru[[ii]]$properties["area"] <- 0
                hru[[ii]]$properties["width"] <- as.numeric( private$shp$width[ii] )
                hru[[ii]]$properties["Dx"] <- as.numeric( private$shp$length[ii] )
                hru[[ii]]$properties["gradient"] <- as.numeric( private$shp$slope[ii] )
                hru[[ii]]$class <- as.list( private$shp[ii,class_names] )
            }

            ## loop the channel routing
            for(ii in 1:nrow(channel_routing)){
                jj <- channel_routing[ii,1]
                kk <- paste(channel_routing[ii,2])
                ff <- channel_routing[ii,3]
                hru[[jj]]$sf_flow_direction[kk] <- ff ## since everything appears once in routing
            }
            ## loop to get outlets
            outlets <- list()
            cnt <- 1
            for(ii in 1:n_channel){
                if( length(hru[[ii]]$sf_flow_direction)==0 ){
                    outlets[[cnt]] <- data.frame(name = paste0("q_sf_",ii),
                                                 id = as.integer( ii ),
                                                 flux = "q_sf", scale = 1.0)
                }
            }
            outlets <- do.call(rbind,outlets)

            ## loop the channel cells
            chn_idx <- which( is.finite( M[["channel"]] ) )
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
            ## pass one to do some properties and rainfall
            for(rw in 1:nrow(flow_routing)){
                ## unpack flow routing
                ii <- flow_routing[rw,1] ## cell number
                gcl <- flow_routing[rw,2]
                dcl <- flow_routing[rw,3]
                w <- flow_routing[rw,4:11]

                jj <- hru_map[ii] ## current hru
                hru[[jj]]$id <- jj
                hru[[jj]]$class <- as.list(json[json[[class_lyr]]== M[[class_lyr]][ii],,drop=FALSE])

                
                ## work out flow routing
                is_ds <- w>0
                if( !any( is_ds ) ){ stop("Hillslope cell with no outflow!") }
                kk <- (ii + delta)[ is_ds ] ## downstream cells
                stopifnot(
                    "All hillslope cells must flow to those with lower id's" = all( jj>hru_map[kk] ),
                    "All hillslope cells must flow to those with bands" = all( M[["band"]][ii]>M[["band"]][kk] )
                )                
                ngh <- paste( hru_map[ kk ] )
                new_ngh <- setdiff( ngh, names( hru[[jj]]$sf_flow_direction ) )
                hru[[jj]]$sf_flow_direction[ new_ngh ] <- 0
                hru[[jj]]$sf_flow_direction[ ngh ] <- hru[[jj]]$sf_flow_direction[ ngh ] + gcl*w[is_ds]

                ## process the other properties
                hru[[jj]]$properties["width"] <- hru[[jj]]$properties["width"] + dcl
                hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] + 1
                hru[[jj]]$properties["gradient"] <- hru[[jj]]$properties["gradient"] + M[["gradient"]][ii]
                hru[[jj]]$band <- min( hru[[jj]]$band , M[["band"]][ii] )

                ## process the rain and pet data
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
            ## pass over hillslope HRU to tidy up properties
            for(jj in (n_channel+1):n_hru){
                hru[[jj]]$properties["gradient"] <- hru[[jj]]$properties["gradient"] / hru[[jj]]$properties["area"]
                hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] * cell_area
                hru[[jj]]$properties["Dx"] <- hru[[jj]]$properties["area"] / hru[[jj]]$properties["width"]
                if( length(hru[[jj]]$sf_flow_direction)==0 ){
                    stop("Hillslope HRU with no outflow")
                }
                hru[[jj]]$sf_flow_direction <- hru[[jj]]$sf_flow_direction / sum(hru[[jj]]$sf_flow_direction)
                kk <- hru[[jj]]$class
            }
                
            ## ############################################
            if(verbose){ cat("Finalising HRUs","\n") }
            for(ii in 1:n_hru){
                hru[[ii]]$id <- as.integer(hru[[ii]]$id - 1) ## since 0 indexed in dynatop
                hru[[ii]]$precip <- list(name = paste0(rainfall_label,names(hru[[ii]]$precip)),
                                         fraction = as.numeric( hru[[ii]]$precip / sum(hru[[ii]]$precip) ))
                hru[[ii]]$pet <- list(name = paste0(pet_label,names(hru[[ii]]$pet)),
                                      fraction = as.numeric( hru[[ii]]$pet / sum(hru[[ii]]$pet) ))
                hru[[ii]]$sf_flow_direction <- list(
                    id = as.integer( as.integer(names(hru[[ii]]$sf_flow_direction)) - 1 ), # since 0 indexed in dynatop
                    fraction = as.numeric( hru[[ii]]$sf_flow_direction ))
                hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
            }
            hru_map[] <- as.integer(hru_map-1)
            outlets$id <- as.integer( outlets$id - 1 )
            outlets$name <- paste0("q_sf_",outlets$id)
            
            ## make output
            if(verbose){ cat("Making output","\n") }
            rst <- terra::rast( private$brk[[ "channel" ]],
                               names = "hru", vals=hru_map)
            model <- list(hru=hru, output_flux = outlets)
            model$map <- paste0(layer_name,".tif")
            terra::writeRaster(rst,model$map,overwrite=TRUE)
            saveRDS(model,paste0(layer_name,".rds"))
        }
    )
)

