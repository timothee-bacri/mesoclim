% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatialdownscale.R
\name{relhumdownscale}
\alias{relhumdownscale}
\title{Downscale relative humidity accounting for downscaled temperature}
\usage{
relhumdownscale(rh, tcc, tcf, dtmc, rhmin = 0)
}
\arguments{
\item{rh}{a coarse-resolution array of relative humidities (percentage)}

\item{tcc}{a coarse-resolution array of temperatures (deg C)}

\item{tcf}{a fine-resolution SpatRast of temperatures (deg C) as returned by
\code{\link[=tempdownscale]{tempdownscale()}}}

\item{dtmc}{a coarse-resolution SpatRast of elevations matching
the resolution, extent and coordinate reference system of \code{rh} and \code{tc}.}

\item{rhmin}{optional single numeric value indicating the minimum realistically
attainable relative humidity within the study region (percentage).}
}
\value{
a stacked SpatRast of relative humidities (percentage) matching the
resolution, coordinate reference system and extent of \code{tcf}.
}
\description{
The function \code{relhumdownscale} is used to spatially downscale relative
}
\examples{
climdata<- read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
dtmf<-terra::rast(system.file('extdata/dtms/dtmf.tif',package='mesoclim'))
tcf<-tempdownscale(climdata,sst=NA,dtmf=dtmf,cad=FALSE,coastal=FALSE)
rhf<-relhumdownscale(climdata$relhum,climdata$temp,tcf,climdata$dtm)
plot_q_layers(rhf)
}
\seealso{
\code{\link[=tempdownscale]{tempdownscale()}}
}
\keyword{spatial}
