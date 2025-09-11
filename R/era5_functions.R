#' @title convert era4 ncdf4 file to format required for model
#' @description The function `era5toclimarray` converts data in a netCDF4 file returned
#' by mcera5 pkg function request_era5 to the correct formal required for subsequent modelling.
#'
#' @param ncfile character vector containing the path and filename of the nc file
#' @param dtmc a SpatRaster object of ERA5 elevations of same or larger extent as aoi and same resolution as ncfile (see details)
#' @param lsm  a SpatRaster object of ERA5 landsea mask of same or larger extent as aoi and same resolution as ncfile
#' @param aoi a  SpatRaster, sf or vect of the area to which outputs are cropped & reprojected - if NA extent/projection same as ncfile
#' @param dtr_cor_fac numeric value to be used in the diurnal temperature range
#' correction of coastal grid cells. Default = 1.285, based on calibration against UK Met Office
#' observations. If set to zero, no correction is applied.
#' @param toArrays logical determining if climate data returned as list of arrays. If FALSE returns list of Spatrasters.
#' @param zo= height of wind speed above ground in metres to be output
#' @return a list of the following:
#' \describe{
#'    \item{dtmc}{Digital elevation of downscaled area in metres (as Spatraster)}
#'    \item{tme}{POSIXlt object of times corresponding to climate observations}
#'    \item{windheight_m}{Height of windspeed data in metres above ground (as numeric)}
#'    \item{tempheight_m}{Height of temperature data in metres above ground (as numeric)}
#'    \item{temp}{Temperature (deg C)}
#'    \item{relhum}{Relative humidity (Percentage)}
#'    \item{pres}{Sea-level atmospheric pressure (kPa)}
#'    \item{swrad}{Total downward shortwave radiation (W/m^2)}
#'    \item{difrad}{Downward diffuse radiation (W / m^2)}
#'    \item{lwrad}{Total downward longwave radiation (W/m^2)}
#'    \item{windspeed}{at 2m above ground (m/s)}
#'    \item{winddir}{Wind direction (decimal degrees)}
#'    \item{prec}{Precipitation (mm)}
#'  }
#' @export
#' @keywords preprocess era5
#' @details the model requires that input climate data are projected using a coordinate reference
#' system in which x and y are in metres. The output data are reprojected using the coordinate reference system and
#' extent of aoi (but SHOULD retain the approximate original grid resolution of the input climate data).
#' Suitable lsm can be derived from ERA5 landsea mask and dtmc can be derived from ERA5 geopotential ancillary data.
#' @examples
#'  \dontrun{
#' ncfile<-'path_to_downloaded_era5_file'
#' aoi<-terra::vect(terra::ext(-7.125,-2.875,49.375,51.625),crs='EPSG:4326')
#' dtm<-terra::rast(system.file('extdata/dtms/era5dtm.tif',package='mesoclim'))
#' era5input<-era5toclimarray(ncfile, dtm=NA, aoi=aoi)
#' plot_q_layers(.rast(era5input$temp,era5input$dtm))
#' checkinputs(era5input,'hour')
#' }
era5toclimarray <- function(ncfile, dtmc, lsm, aoi=NA, dtr_cor_fac = 1.285, toArrays=TRUE, zo=10)  {
  # Get extent of aoi and check equal or smaller than dtmc and lsm inputs
  if (class(aoi)[1] != "logical"){
    if (!class(aoi)[1] %in% c("SpatRaster", "SpatVector",
                              "sf"))
      stop("Parameter aoi NOT of suitable spatial class ")
    if (class(aoi)[1] == "sf") aoi <- vect(aoi)
    # Check dtmc and lsm extent equals or exceeds aoi
    if(ext(dtmc)<ext(project(aoi,crs(dtmc)))) stop("dtmc smaller than aoi!!!")
    if(ext(lsm)<ext(project(aoi,crs(dtmc)))) stop("dtmc smaller than aoi!!!")
  }
  # If no aoi assume same as dtmc
  if(class(aoi)[1]=="logical") aoi <- dtmc
  units(dtmc)<-'m'
  names(dtmc)<-'Elevation'

  # Check vars and get time
  nc <- nc_open(ncfile)
  era5vars <- names(nc$var)
  if("expver" %in% era5vars) tme<-as.POSIXlt(nc$var$expver$dim[[1]]$vals,tz='GMT',origin="1970-01-01")
  if(!"expver" %in% era5vars) tme<-as.POSIXlt(time(rast(ncfile,"t2m")), tz = "UTC")
  nc_close(nc)


  t2m <- rast(ncfile, subds = "t2m") %>% crop(dtmc)
  d2m <- rast(ncfile, subds = "d2m") %>% crop(dtmc)
  if('sp' %in% era5vars) pres <- rast(ncfile, subds = "sp") %>% crop(dtmc) else if('msl' %in% era5vars){
    psl <- rast(ncfile, subds = "msl")
    pres<-psl * (((293-0.0065*dtmc)/293)^5.26) # convert to sea level pressure
  } else stop("Missing necessary pressure variables!!!")

  u10 <- rast(ncfile, subds = "u10") %>% crop(dtmc)
  v10 <- rast(ncfile, subds = "v10") %>% crop(dtmc)
  tp <- rast(ncfile, subds = "tp") %>% crop(dtmc)
  msdwlwrf <- rast(ncfile, subds = "msdwlwrf") %>% crop(dtmc)
  if (all(c("fdir", "ssrd") %in% era5vars)) {
    fdir <- crop(rast(ncfile, subds = "fdir")/3600,dtmc)
    ssrd <- crop(rast(ncfile, subds = "ssrd")/3600,dtmc)
  } else if (all(c("msdwswrf", "msdrswrf") %in% era5vars)) {
    fdir <- rast(ncfile, subds = "msdwswrf") %>% crop(dtmc)
    ssrd <- rast(ncfile, subds = "msdrswrf") %>% crop(dtmc)
  } else stop("Missing necessary SW radiation variables!!!")

  # Aggregate dtmc if not of same resolution as ncfile data
  if(res(dtmc)[1]!=res(t2m)[1]){
    agf <- terra::res(t2m)[1]/terra::res(dtmc)[1]
    dtmc <- terra::aggregate(dtmc, fact = agf, fun = mean, na.rm = T)
  }

  # Convert temp to Celsius
  t2m <- t2m - 273.15

  # Coastal correction of temperature
  lsm_e<-crop(lsm,dtmc)
  tmn <- .ehr(.hourtoday(as.array(t2m), min))
  a<-as.array(rep(lsm_e,dim(tmn)[3]))
  mu<-(1-a)*dtr_cor_fac+1
  tc <- .rast(((as.array(t2m)) - tmn) * mu + tmn, t2m)

  # Reproject and crop all variables to aoi
  tc <- terra::project(tc, crs(aoi))
  d2m<- terra::project(d2m,crs(aoi))
  pres <- terra::project(pres, crs(aoi))
  u10 <- terra::project(u10, crs(aoi))
  v10 <- terra::project(v10, crs(aoi))
  tp <- terra::project(tp, crs(aoi))
  msdwlwrf <- terra::project(msdwlwrf, crs(aoi))
  fdir <- terra::project(fdir, crs(aoi))
  ssrd <- terra::project(ssrd, crs(aoi))
  dtmc<-project(dtmc,crs(aoi))
  lsm<-project(lsm,crs(aoi))

  # Calculate output variables
  ea <- .rast(.satvap(as.array(d2m)-273.15), tc)
  temp <- as.array(tc)
  relhum <- (as.array(ea)/.satvap(temp)) * 100
  pres <- as.array(pres)/1000 ## Surface Pressure!!!!
  swrad <- as.array(ssrd)
  difrad <- swrad - as.array(fdir)
  lwrad <- as.array(msdwlwrf)
  windspeed <- sqrt(as.array(u10)^2 + as.array(v10)^2) * log(67.8 * zo - 5.42)/log(67.8 * 10 - 5.42)
  winddir <- as.array((terra::atan2(u10, v10) * 180/pi + 180)%%360)
  prec <- as.array(tp) * 1000

  # Format outputs
  out <- list(dtm = dtmc, tme = tme, windheight_m = zo, tempheight_m = 2)
  climout <- list(temp = temp, relhum = relhum, pres = pres,
                  swrad = swrad, difrad = difrad, lwrad = lwrad, windspeed = windspeed,
                  winddir = winddir, prec = prec)
  climunits<-c('degC','%','kPa','watt/m^2','watt/m^2','watt/m^2','m/s','deg','mm')
  # Round values
  climround<-c(2,1,1,1,1,1,3,1,4)
  for(n in 1:length(climout)){
    v <- names(climout)[n]
    climout[[v]]<-round(climout[[v]],climround[n])
  }
  # For spatrasters add units and time
  if (toArrays == FALSE){
    climout <- lapply(climout, .rast, tem = dtmc)
    for(n in 1:length(climout)){
      v <- names(climout)[n]
      u<-climunits[n]
      terra::time(climout[[v]])<-tme
      names(climout[[v]])<-tme
      units(climout[[v]])<-u
    }
  }
  out <- c(out, climout)
  return(out)
}
