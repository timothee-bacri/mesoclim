library(parallel)
detectCores()

t0<-now()
mesoclimate <- spatialdownscale(
  climdata,
  sstdata,
  dtmf, dtmm,
  # NEW PAAMETERS

  cad = TRUE, coastal = TRUE, thgto = 2, whgto =   2,
  # NEW PARAMETER - outputs tmean option
  include_tmean=TRUE,
  rhmin = 20, pksealevel = TRUE, patchsim = FALSE,
  terrainshade = FALSE, precipmethod = "Elev", fast = TRUE, noraincut = 0.01
)
tend<-now()-t0
print(tend) # 37 secs for 1 year


clim.list<-list()
for (n in 1:12){
  print(n)
  smonth<-climdata$tme[which(month(climdata$tme)==n)][[1]]
  emonth<-rev(climdata$tme[which(month(climdata$tme)==n)])[[1]]
  clim.list[[n]]<-subset_climdata(climdata,smonth,emonth)
}

wca<-calculate_windcoeffs(dtmc,dtmm,dtmf,zo=2)

# Cold air drainage basins - as above and ONLY if using coastal correction - can take several minutes for large areas
basins<-basindelin(dtmf, boundary = 2)

#terrain shading
results<-calculate_terrain_shading(dtmf,steps=24,toArrays=FALSE)
skyview<-results$skyview
horizon<-results$horizon

t0<-now()
meso.list<-lapply(clim.list,spatialdownscale,
                    sst=sstdata, dtmf=dtmf, dtmm =dtmm,
                  basins=basins, wca=wca, skyview=skyview,horizon=horizon,
                  cad = TRUE,coastal = TRUE, thgto =2, whgto=2,
                    include_tmean=TRUE,rhmin = 20, pksealevel = TRUE,
                    patchsim = FALSE, terrainshade = TRUE,
                    precipmethod = "Elev",fast = TRUE, noraincut = 0.01)
tend<-now()-t0
print(tend)  # 3.8 mins when preprocessing topog


t0<-now()
meso.list<-mclapply(clim.list[1:8],spatialdownscale,
                    sst=sstdata, dtmf=dtmf, dtmm =dtmm,
                    cad = TRUE,coastal = TRUE, thgto =2, whgto=2,
                    include_tmean=TRUE,rhmin = 20, pksealevel = TRUE,
                    patchsim = FALSE, terrainshade = TRUE,
                    precipmethod = "Elev",fast = TRUE, noraincut = 0.01,
                    mc.preschedule = FALSE)
tend<-now()-t0
print(tend) # 3.106675 mins 8 months = 2.75mins
res<-mccollect(meso.list)
