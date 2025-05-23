% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/workerfunctions.R
\name{.sea_to_coast}
\alias{.sea_to_coast}
\title{Interpolates sea to coastal sea surface temperature
Resamples sst data to same extent, resolution and projection as aoi
Interpolates from mean of adjacent sea cells}
\usage{
.sea_to_coast(sst.r, aoi.r, ext_cells = 8)
}
\arguments{
\item{sst.r}{\itemize{
\item sea surface temperature data
}}

\item{aoi.r}{\itemize{
\item dtm of area of interest where sea = NA
}}
}
\value{
sst spatrasters of same projection, res and extent as aoi with matching coastline
}
\description{
Interpolates sea to coastal sea surface temperature
Resamples sst data to same extent, resolution and projection as aoi
Interpolates from mean of adjacent sea cells
}
\examples{
sst.r<-rast(ukcp18sst[[1]])
aoi.r<-
}
