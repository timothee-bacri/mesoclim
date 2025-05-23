% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_functions.R
\name{read_climdata}
\alias{read_climdata}
\title{Read climate data}
\usage{
read_climdata(filepath)
}
\arguments{
\item{filepath}{\itemize{
\item either a pre-loaded list or filepath to a list of prepared climate data as written by \code{write_climdata()}
}}
}
\value{
R list of climate data
}
\description{
Reads a climate dataset as output by for example \code{ukcp18toclimdata()} and written by \code{write_climdata()}
Unwraps any PackedSpatRasters to Spatrasters
}
\examples{
climdata<-read_climdata(system.file('extdata/preprepdata/ukcp18rcm.Rds',package='mesoclim'))
climdata<-read_climdata(ukcpinput)
}
\keyword{data}
\keyword{preprocess}
