\name{check.temporal}
\alias{check.temporal}
\title{Check the temporal coordinate of a point pattern}
\description{
  A function to check the temporal coordinate of a spatio-temporal 
  point pattern to be in ascending order. 
}
\usage{
  check.temporal(X)
}
\arguments{
  \item{X}{Marked spatio-temporal point pattern (object of class
    \code{"ppx"}).}
}
\value{
  A character string.
}
\details{
  The function produces a message (character string) indicating 
  whether the temporal coordinate of the point pattern is in 
  asceding order or not.
}
\seealso{
  \code{\link{ppx}}.
}
\examples{
  check.temporal(iran.quakes)

  check.temporal(jap.quakes)
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
  \url{http://www.razi.ac.ir/jalilian/}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
