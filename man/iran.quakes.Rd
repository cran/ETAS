\name{iran.quakes}
\alias{iran.quakes}
\docType{data}
\title{20th century earthquakes of Iran}
\description{
  A point pattern giving occurrence time (in days), location of epicenter
  (longitude and latitude both in decimal degrees) and magnitude of
  shalow earthquakes (depth < 100 km) occurred since 1902-01-01 till
  2012-06-01 in Iran and its vicinity. Earthquakes with magnitude less
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
      \code{mag.type} \tab Magnitude type \cr
      \code{depth} \tab Depth (km)\cr
      \code{ref}   \tab Reference \cr
      \code{date}  \tab Actual date and time (UTC)\cr
      }
  \item and domian of data
      \tabular{ll}{
      \code{time domain} \tab [0, 15492] : from 1970-01-01 (time origin 0)
      to 2012-06-01 (15492) \cr
      \code{spatial domain} \tab [41, 67] * [20.5, 44.5] \cr
      \code{mark domain}  \tab [4.5, \code{Inf}) \cr
      }
  }
}
\usage{data(iran.quakes)}
\source{International Institute of Earthquake Engineering and Seismology (IIEES): www.iiees.ac.ir/}
\section{Notes}{
more points about the data, their accuracy and ....
}
\references{
  Jalilian, A. (2012). Modeling Earthquakes of Iran.
}
\keyword{datasets}
\keyword{spatial}
\keyword{earthquak}
