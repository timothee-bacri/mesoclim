#' @title Downscale longwave radiation
#' @description Downscales an array of longwave downward radiation data with option to
#' simulate terrain shading of skyview
#'
#' @param lwrad - coarse resolution downward longwave radiation (W/m^2)  as spatraster or array.
#' @param tc - coarse resolution temperature (deg C) as spatraster or array.
#' @param tcf - fine resolution spatraster of temperature (deg C) down-scaled to resolution of `dtmf`.
#' @param tme POSIXlt object of times corresponding to radiation values in
#' `lwrad`, `tc` and `tcf`.
#' @param dtmf a fine-resolution SpatRast of elevations.
#' @param dtmc SpatRast of elevations matching the resolution, extend and coordinate
#' reference system of `lwrad` and `tc`.
#' @param terrainshade - boolean TRUE/FALSE
#'
#' @return  a multi-layer SpatRast of downward longwave radiation (W/m^2)
#' matching the resolution of dtmf.
#' @details radiation is downscaled by computing from temperature the difference in upward
#' longwave radiation at fine and coarse resolutions. Corrected LW down is then
#' corrected for terrain shading of skyview if requested.
#' @export
#' @examples
#' climdata<- read_climdata(mesoclim::ukcpinput)
#' dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
#' dtmm<-terra::rast(system.file('extdata/dtms/dtmm.tif',package='mesoclim'))
#' tmean_c<-(climdata$tmin+climdata$tmax)/2
#' wsf<-winddownscale(climdata$windspeed, climdata$winddir, dtmf, dtmm, climdata$dtm, zi=climdata$windheight_m)
#' dailytemps<-tempdaily_downscale(climdata,NA,terra::unwrap(mesoclim::ukcp18sst),dtmf,dtmm,NA,wsf,cad = FALSE,coastal = FALSE,2,2)
#' lw<-lwdownscale(climdata$lwrad,tmean_c,dailytemps$tmean, climdata$tme, dtmf=dtmf, dtmc=climdata$dtm)
#' terra::panel(c(lw[[1]],lw[[5]],lw[[10]],lw[[20]]),main=paste0("LW down ",c(1,5,10,15),"/05/2018"))
#'
lwdownscale_original<-function(lwrad, tc, tcf, tme, dtmf, dtmc, skyview=NA, terrainshade = TRUE) {
  if(class(lwrad)[1]=="array") lwrad<-.rast(lwrad,dtmc)
  lwf<-.resample(lwrad,dtmf, msk=TRUE)
  lwupc<-.lwup(tc)
  if(class(lwupc)[1]=="array") lwupc<-.rast(lwupc,dtmc)
  lwupc<-.resample(lwupc,dtmf)
  lwupf<-.lwup(tcf)
  lwf<-lwf+(lwupc-lwupf)
  if (terrainshade) {
    if(inherits(skyview,'logical')) skyview<-.skyview(dtmf)
    #lwf<-.rast(.is(lwf)*skyview,dtmf)
    lwf<-lwf*skyview
  }
  return(lwf)
}

lwdownscale_new<-function(lwrad,  method=c("interp","idso","dob"),
                          tcf, rhf, pkf, tme, dtmf, dtmc,
                          skyview=NA, terrainshade = TRUE,
                          cloudcover=NA, cloudcorrect=FALSE) {
  if(method=="interp"){
    if(class(lwrad)[1]=="array") lwrad<-.rast(lwrad,dtmc)
    lwf<-.resample(lwrad,dtmf, msk=TRUE)
  }

  if(method=="idso"){
    sb<-5.67e-08
    ea<-converthumidity(rhf, tc=tcf,pk=pkf)
    em<-0.7 + 5.95e-04 * ea * exp(1500/(tcf+273.15))
    lwf<-em*sb*(tcf+273.15)^4
  }

  if(method=="dob"){
    Tref<-273.15
    wref<-25
    ea<-converthumidity(rhf, tc=tcf,pk=pkf)
    w = 465 * (ea / tcf)
    lwf<-59.38 + 113.7 * ((tcf+273.15) / Tref)^6 + 96.96 * sqrt(w/wref)
  }

  if(cloudcorrect){
    # add cloud cover correction if not provided
    sb<-5.67e-08
    cld_hght<-1000
    ec<-1
    ea<-converthumidity(rhf, tc=tcf,pk=pkf)
    em<-0.7 + 5.95e-04 * ea * exp(1500/(tcf+273.15))
    Tc<-(273.15+tcf) - 0.0065*(cld_hght-dtmf)
    f8c<- -0.6732 + 0.006240 * Tc - 0.9140e-05 * Tc^2
    e8z<- 0.24 +2.98e-06 * em^2 * exp(3000/(tcf+273.15))
    e8<-e8z*(1.4-0.4*e8z)
    t8<-1-e8
    lwfc<-t8 * cloudcover * ec * f8c * sb * Tc^4
    lwf<-lwf+lwfc
  }

  if (terrainshade) {
    if(inherits(skyview,'logical')) skyview<-.skyview(dtmf)
    #lwf<-.rast(.is(lwf)*skyview,dtmf)
    lwf<-lwf*skyview
  }

  return(lwf)
}

lwfinterp<-lwdownscale_new(climdata$lwrad,  method="interp",
                           tcf=mesoclimate$tmean, rhf=mesoclimate$relhum,
                           pkf=mesoclimate$pres, tme=mesoclimate$tme,
                           dtmf=mesoclimate$dtm, dtmc=climdata$dtm,
                           skyview=NA, terrainshade = FALSE,
                           cloudcover=NA, cloudcorrect=FALSE)
lwfidso<-lwdownscale_new(climdata$lwrad,  method="idso",
                           tcf=mesoclimate$tmean, rhf=mesoclimate$relhum,
                           pkf=mesoclimate$pres, tme=mesoclimate$tme,
                           dtmf=mesoclimate$dtm, dtmc=climdata$dtm,
                           skyview=NA, terrainshade = FALSE,
                           cloudcover=NA, cloudcorrect=FALSE)
lwfdob<-lwdownscale_new(climdata$lwrad,  method="dob",
                           tcf=mesoclimate$tmean, rhf=mesoclimate$relhum,
                           pkf=mesoclimate$pres, tme=mesoclimate$tme,
                           dtmf=mesoclimate$dtm, dtmc=climdata$dtm,
                           skyview=NA, terrainshade = FALSE,
                           cloudcover=NA, cloudcorrect=FALSE)

r<-c(lwfinterp[[1]],lwfidso[[1]],lwfdob[[1]])
names(r)<-c("interp","idso","dob")
plot(r)

cloudcover<-resample(.rast(climdata$cloud,climdata$dtm),dtmf)/100
lwfinttc<-lwdownscale_new(climdata$lwrad,  method="interp",
                           tcf=mesoclimate$tmean, rhf=mesoclimate$relhum,
                           pkf=mesoclimate$pres, tme=mesoclimate$tme,
                           dtmf=mesoclimate$dtm, dtmc=climdata$dtm,
                           skyview=NA, terrainshade = TRUE,
                           cloudcover=NA, cloudcorrect=FALSE)
lwfidcc<-lwdownscale_new(climdata$lwrad,  method="idso",
                         tcf=mesoclimate$tmean, rhf=mesoclimate$relhum,
                         pkf=mesoclimate$pres, tme=mesoclimate$tme,
                         dtmf=mesoclimate$dtm, dtmc=climdata$dtm,
                         skyview=NA, terrainshade = FALSE,
                         cloudcover=cloudcover, cloudcorrect=TRUE)
lwfdobcc<-lwdownscale_new(climdata$lwrad,  method="dob",
                        tcf=mesoclimate$tmean, rhf=mesoclimate$relhum,
                        pkf=mesoclimate$pres, tme=mesoclimate$tme,
                        dtmf=mesoclimate$dtm, dtmc=climdata$dtm,
                        skyview=NA, terrainshade = FALSE,
                        cloudcover=cloudcover, cloudcorrect=TRUE)

r2<-c(lwfinttc[[1]],lwfidcc[[1]],lwfdobcc[[1]])
names(r2)<-c("interp_tc","idso_cc","dob_cc")
plot(r2)
