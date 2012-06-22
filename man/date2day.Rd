\name{date2day}
\alias{date2day}
\title{Convert date-time data to numeric data in days}
\description{
  A function to convert date-time data to days with respect to a date-time origin. 
}
\usage{
  date2day(dates, start, tz="", \dots)
}
\arguments{
  \item{dates}{A date-time or date object. Typically, it is a character 
    vector containing date-time information.}
  \item{start}{A date-time or date object. Determines the origin of the conversion.}
  \item{tz}{Optional. Timezone specification to be used for the conversion.}
  \item{\dots}{Arguments to be passed to \code{as.POSIXlt}.}
}
\value{
  A numeric vector of the same length as \code{dates}.
}
\details{
  The arguments \code{dates} and \code{start} must be of 
  appropriate format to be passed to \code{as.POSIXlt} function.
}
\seealso{

  \code{\link{as.POSIXlt}} and \code{\link{difftime}} for appropriate format of the data
  to be converted.
}
\examples{
  # date-time data of Iran's earthquakes between 1902/0/01 and 2012/06/01
  dt <- iran.quakes$data$date
  # origin of the conversion
  start <- "1902/01/01 00:00:00"
  # time in days since 1902-01-01 (UTC)
  date2day(dt, start, tz="GMT")
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
  \url{http://www.razi.ac.ir/jalilian/}
}
\keyword{spatial}
\keyword{math}
\keyword{date time}
