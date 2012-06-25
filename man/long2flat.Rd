\name{long2flat}
\alias{long2flat}
\title{Project longitude-latitude coordinates to flat map coordinates}
\description{
  A function to adjust the longitude of coordinates of a set of points 
  in a geographical region on earth in order to have a flat map coordinates. 
}
\usage{
   long2flat(X, win)
}
\arguments{
  \item{X}{Marked spatio-temporal point pattern (object of class
    \code{ppx}).}
  \item{win}{Target geographical region of study (object of class
    \code{owin}).}
}
\value{
  A list with components
  \describe{
    \item{X:}{The same point pattern as the input argument \code{X} except that
    the longitude coordinates of its points are adjusted.}
    \item{win:}{The same window as the input argument \code{win} except 
    that the longitude coordinates of its vertexes are adjusted.}
  }
}
\details{
  The algorithm does not change the latitude of points but adjust the longitude 
  of points with respect to the centroid of the target geographical region \code{win}
  in order to get coordinates on a flat surface. 
}
\seealso{

  \code{\link{ppx}}

  \code{\link{owin}}
}
\examples{
  # specifiying the geographical region
  win <- owin(c(41, 66), c(24, 42))
  \dontrun{
  plot(iran.quakes$data[, 2:3])
  plot(win, add=TRUE)
  }
  # projecting log-lat coordinates into flat map coordinates
  proj <- long2flat(iran.quakes, win)
  \dontrun{
  plot(proj$X$data[, 2:3])
  plot(proj$win, add=TRUE)
  }
   
  # specifiying the geographical region
  jwin <- owin(poly=list(x=c(134.0, 137.9, 143.1, 144.9, 
                        147.8, 137.8, 137.4, 135.1, 130.6),
                         y=c(31.9, 33.0, 33.2, 35.2, 41.3, 
                             44.2, 40.2, 38.0, 35.4)))
  \dontrun{
  plot(jap.quakes$data[, 2:3])
  plot(jwin, add=TRUE)
  }
  # projecting log-lat coordinates into flat map coordinates
  proj <- long2flat(jap.quakes, jwin)
  # specifying time period
  \dontrun{
  plot(proj$X$data[, 2:3])
  plot(proj$win, add=TRUE)
  }
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
  \url{http://www.razi.ac.ir/jalilian/}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
