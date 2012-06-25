\name{jap.quakes}
\alias{jap.quakes}
\docType{data}
\title{Earthquakes of Japan}
\description{
  A point pattern giving occurrence time (in days), location of epicenter
  (longitude and latitude both in decimal degrees) and magnitude of
  shalow earthquakes (depth < 100 km) occurred since 1926-01-08 till
  2008-01-01 in Japan and its vicinity. Earthquakes with magnitude less
  than 4.5 have been excluded from data.
} 
\format{
  An object of class \code{"ppx"} representing informations about 
  \itemize{
  \item data
    \tabular{ll}{
      \code{time}  \tab Occurrence time (days)\cr
      \code{long}  \tab Longitude of epicenter (degree)\cr
      \code{lat}   \tab Latitude of epicenter (degree)\cr
      \code{mag}   \tab Magnitude (different scales)\cr
      \code{depth} \tab Depth (km)\cr
      \code{date}  \tab Actual date  (local time) \cr
      }
  \item and domian of data
      \tabular{ll}{
      \code{time domain}    \tab [0, 15492] : from 1926-01-08 (time origin 0)
      to 2007-12-29 (29943) \cr
      \code{spatial domain} \tab [128, 145] x [27, 45] \cr
      \code{mark domain}  \tab [4.5, Inf) \cr
      }
  }
}
\usage{data(jap.quakes)}
\source{(?)}
\section{Notes}{
more points about the data, their accuracy and ....
}
\references{
 Zhuang and Ogata (?)
} 
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
