#!/usr/bin/env Rscript

## Run this as Rscript
## Rscript --vanilla RunForecast.r '2022-01'

### --- Load packages
library(optparse, quietly = TRUE)

## Create command line options
option_list <- list(
    make_option(c("-d", "--date"), action = "store", default = "2022-01",
                help = "Date specified as string YYYY-mm. Default is 2022-01"),
    make_option(c("-m", "--model"), action = "store",
                default = "v2-tc2h_nb10_res5_captures",
                help = "Name of model to use"),
    make_option(c("-o", "--out"), action = "store",
                default = "Predictions",
                help = "Path where forecasts are saved. Default is /Predictions"),
    make_option(c("-v", "--verbose"), action = "store",
                default = TRUE,
                help = "Prints detailed output")
)

## Parse and store parameters
opts <- parse_args(OptionParser(option_list=option_list))


### *****************************

## Uses fitted models and monthly precipitation averages to forecast
## rodent trap success or viral prevalence three months ahead of the
## most recent rainfall data. 

## Outline of code
## Step 1: Read in the raster data
## Step 2: Apply model to weather time-series to calculate trap success
## Step 3: Multiple trap success by the LASV pathogen layer to calculate LASV spillover risk

### *****************************

verbose <- opts$verbose

if(verbose){writeLines('\n \n *** Loading packages *** ')}
library(lubridate, quietly = TRUE)
library(xgboost, quietly = TRUE)
library(chirps, quietly = TRUE)
library(terra, quietly = TRUE)

## Get forecast date, check that it is valid



forecast.date <- opts$date
forecast.date = paste0(forecast.date, '-01')
forecast.date <- as.Date(forecast.date)
if(!is.Date(forecast.date)){
    writeLines("Provide date in YYYY-mm format; e.g., '2020-01'")
}

### --- Default settings

prefix <-  opts$model ## Name of model (prefix set in Main_Fit_XGB.r)
include.types <- 1 ## Include Type == 'House' only

## Specify where the forecast rasters will be saved
full.forecast.path <- opts$out 


## - Specify folder that contains predictor stacks
raster.data.path <- "Data/Rasters"

## - Specify location of MODIS monthly precipitation rasters
prec.path <- 'Data/Rasters/Precipitation'##opts$prec

## - Specify shapefile (West Africa) path and layer name
shapefile.path <- paste('Data/Shapefiles/West_Africa', sep = '')
shapefile.name <- 'foc'

### --- END USER INPUT

if(verbose){writeLines('\n \n *** Reading in shapefiles *** ')}

## Load in shapefile of west african countries
foc.shp <- vect(paste0(shapefile.path, '/', shapefile.name, '.shp'))

## - Location of mastomys natalensis rangemap for masking
masto.rangemap <- vect('Data/Shapefiles/Mn_Rangemap/data_0.shp')

## Define locations of model data and models
models.info.path <-  paste0('Models/Rodent_Layer/', prefix, '/')
models.path <-paste0(models.info.path, 'Model_Fits')
data.path <- paste0(models.info.path, 'Data_Files')


## --- Get model names that must be loaded in
if(verbose){writeLines('\n \n *** Reading in model parameters *** ')}

## Load in full hyperparameter set that was run
meta.grid <- read.csv(paste0(models.info.path, 'meta_filled.csv'))
summary.fits <- read.csv(paste0(models.info.path, 'summary_fits.csv'))

jj.meta.grid <- which(meta.grid$test.site=='none')
parms <- unique(meta.grid[jj.meta.grid,])

## --- Delete models that are not required
all.mods <- list.files(paste0("Models/Rodent_Layer/", prefix, '/Model_Fits'), full.names = TRUE)
keep.model.name <- paste0(models.path, '/mod',jj.meta.grid, sep = '')
delete.models <- !(all.mods %in% keep.model.name)
unlink(all.mods[delete.models])

## --- Set up directory system in which to save forecasts

## Define the names of storage directories
raster.path <- paste0(full.forecast.path)
house.tif.path <- paste0(full.forecast.path, '/House/')
house.lasv.tif.path <- paste0(full.forecast.path, '/Forecast/')

## Create directory structure
dir.create(path = house.tif.path, showWarnings = FALSE, recursive = TRUE)
dir.create(path = house.lasv.tif.path, showWarnings = FALSE, recursive = TRUE)

## --- Load in rainfall rasters

if(verbose){writeLines('\n \n *** Reading in raster sets *** ')}

include.prec <- unique(parms$include.prec)
if(length(include.prec) > 1){print('!!!WARNING: Ambiguous parameters'); q()}


if(verbose){writeLines('--Rainfall ')}


dir.create('Data/Rasters/Precipitation/', recursive = TRUE)


prec.stack <- rast()

for(ii in 3:12){
    prec.start.date <- forecast.date - months(ii)
    prec.end.date <- forecast.date - months(ii-1) - days(1)
    date.format <- format(prec.start.date, "%Y.%m")
    
    chirps.month.name = paste0("chirps-v2.0.",
                               format(prec.start.date, "%Y.%m"),
                               ".tif")
    chirps.filename = paste0(prec.path, '/', chirps.month.name)
    
    file.isthere <- file.exists(x = chirps.filename) 

    ## Check if monthly precip datum is present, if not download and save it to file
    if(!file.isthere){
        mean.month.stack <- NULL

        ## AJB: I notice that the Chirps API will not download all of March, 2022, even though
        ## this data is available on the CHC website. The API does allow downloads of
        ## March 1 - March 21, however. Create a while loop that progressively tries date
        ## spans from March 1 - (March 31 - ddates). Ntries and ddates will track how much
        ## of a fudge the resulting monthly average rainfall calculation is. 
        ntries = 0
        ddate <- 0
        success <- FALSE
        while(!success & ntries < 3){
            try({
                ntries = ntries + 1
                out = get_chirps(foc.shp,
                                 dates = c(prec.start.date, prec.end.date - days(ddate)),
                                 server = "CHC",
                                 as.raster = TRUE)
                out[out < 0] = NA
                mean.month.stack <- mean(out, na.rm = TRUE)
            })
            if(is.null(mean.month.stack)){
                ddate = ddate + 5
            }else{
                success = TRUE
            }
        }## End loop trying to acquire rainfall

        ## Let the user know if the rainfall data was not generated. Otherwise,
        ## save the rainfall data ONLY WHEN ddate==0 (meaning that the full
        ## month of rainfall was downloaded). This will force the code to try
        ## again in a month. 
        if(is.null(mean.month.stack)){
            writeLines('*** ERROR: insufficient rainfall data for this prediction')
            q()
        }else{
            if(ddate == 0){
                writeRaster(x = mean.month.stack,
                            file = chirps.filename)
            }
        }
        if(verbose){writeLines(paste0('--- ', date.format, ' processed'))}

    }else{ ## If monthly precip raster is in prec.path then just load it
        mean.month.stack <- rast(chirps.filename)
        if(verbose){writeLines(paste0('--- ', date.format, ' loaded'))}
    }

    prec.stack <- c(prec.stack, mean.month.stack)
        
}## End loop through precipitation lags

prec.names <- paste0('P', 3:12)
names(prec.stack) <- prec.names

## Used to apply a XGB model to a raster
xgb.pred <- function(model, data) {
    predict(model, as.matrix(data))
}


## --- Incorporate environmental and rainfall features into dataset

if(verbose){writeLines('--Habitat ')}
hab.stack <- rast(x = paste0(raster.data.path,
                                       '/Habitat',
                                       "/predictor_stack_05.grd"))
hab.preds <- names(hab.stack)
## Choose extent from habitat as default
ext.obj <- ext(hab.stack)

prec.stack <- resample(prec.stack, hab.stack)
prec.stack <- crop(prec.stack, hab.stack)

## --- Prepare Nights and Type rasters

## Set up "Nights' raster variable
zeros.rast <- hab.stack[[1]]*0.0
zeros.rast <- crop(zeros.rast, ext.obj)
night.rast <- zeros.rast + 1.0 ## Night
names(night.rast) <- 'Night'
night.rast <- resample(night.rast, hab.stack)
night.rast <- crop(night.rast, hab.stack)

## Define seperate feature stacks for house and intown predictions
type.stack <- c(zeros.rast, zeros.rast, zeros.rast)
names(type.stack) <- c('House', 'InTown', 'OutTown')
type.stack <- resample(type.stack, hab.stack)
type.stack <- crop(type.stack, hab.stack)

## Location of LASV raster
if(verbose){writeLines('--Pathogen layer ')}

lasv.location <- 'Models/Pathogen_Layer/Lassa_Layer.tif'
lasv.rast <- rast(lasv.location)
lasv.rast = mask(lasv.rast, masto.rangemap, updatevalue = 0)
lasv.max <- minmax(lasv.rast)[2]
lasv.rast <- resample(lasv.rast, hab.stack)
lasv.rast <- crop(lasv.rast, hab.stack)

## Define dataframe of dates on which to apply the model. A single NA date if model is static. 
grid.dat <- expand.grid(Date = forecast.date,
                        Type = (c('House', 'InTown', 'OutTown')[1:include.types]))

mcl.fun <- function(ti){
    library('terra')
    library('xgboost')
    ## Get focal date of prediction
    cur.date <- grid.dat$Date[ti]
    type <- paste(grid.dat$Type[ti])
    
    ## Join all features together into a single stack
    type.stack[[type]] <- 1.0

    feature.stack <- c(hab.stack, prec.stack, night.rast, type.stack)
    ## Define raster stack variables in which predictions will be saved
    prediction.stack <- rast()
    ## Load models and apply to predictor stack    
    for(mod.i in jj.meta.grid){
        model.name <- paste0(models.path, '/mod',mod.i, sep = '')
        if(file.exists(model.name)){

            ## --- Load model and predict in house and intown sites
            xgb.mod <- xgb.load(model.name)

            pred.names <- unlist(read.table(paste0(data.path,
                                                   '/pred_names_',mod.i, '.txt', sep = '')),
                                 use.names = FALSE)

            pred.names <- pred.names[which((pred.names %in% names(feature.stack)))]      
            rast.out <- predict(feature.stack[[pred.names]], xgb.mod, fun=xgb.pred)
            
            prediction.stack <- c(prediction.stack, rast.out)
        } ## Check if file exists
    }## End loop through boots
    ## Average together raster predictions and calculate standard deviation
    mean.rast <- mean(prediction.stack)

    ti.name.rast <- format(cur.date, '%Y-%m')

    ## Correct water body areas in raster predictions (shows up as NA, should be 0)
    mean.rast[is.na(mean.rast)] <- 0
    mean.rast <- terra::mask(mean.rast, foc.shp, updatevalue = NA)
    ts.minmax <- minmax(mean.rast)
    
    ## Write rasters to file
    tif.path <- ifelse(type=='House', house.tif.path,
                ifelse(type=='InTown', intown.tif.path,
                       outtown.tif.path))
    writeRaster(x = mean.rast,
                filename = paste(tif.path, ti.name.rast, '.tif', sep = ''),
                overwrite = TRUE)    
    
    ## Multiply by LASV layer to calculate spillover risk
    lasv.risk <- mean.rast*lasv.rast
    ## Save a copy lasv.risk raster as tif
    writeRaster(lasv.risk,
                filename = paste0(house.lasv.tif.path, '/', ti.name.rast, '.tif'), 
                overwrite = TRUE)


    
    ## Return min, max of rasters (useful for plotting)
    return(
        data.frame(Date = grid.dat$Date[ti],
                   Type = type,
                   min.val = ts.minmax[1],
                   max.val = ts.minmax[2])
    )
} 

## --- Calculate prediction rasters by calling the above function for each prediction timepoint 

if(verbose){writeLines('\n \n *** Applying models *** ')}

starttime <- Sys.time()
send.vars <- ls(all = TRUE)

out <- lapply(1:nrow(grid.dat), mcl.fun)
if(verbose){writeLines(paste('-Finished. Total Time: ', Sys.time() - starttime))}

##a = rast('Predictions/Forecast/2022-05.tif')

