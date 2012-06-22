\name{print.etas}
\alias{print.etas}
\title{Print a Fitted ETAS Model}
\description{
  Default \code{print} method for a fitted ETAS model.
}
\usage{
 \method{print}{etas}(x,\dots)
}
\arguments{
  \item{x}{
    A fitted ETAS model, obtained from
    the model-fitting algorithm \code{\link{etas}}.
    An object of class \code{"etas"}.
  }
  \item{\dots}{Ignored.}
}
\value{
  None.
}
\details{
  This is the \code{print} method for the class \code{"etas"}.
  It prints information about the fitted model in a sensible format.

}
\seealso{

  \code{\link{etas}} for fitting an ETAS model.
}
\examples{
  win <- owin(c(41, 66), c(24, 42))
  proj <- long2flat(iran.quakes, win)
  param01 <- c(0.4339678,
   0.1988628,
   0.0345206,
   1.6290137,
   1.1286776,
   0.0072539,
   2.1705884,
   0.5706402)
 \dontrun{
  res <- etas(proj$X, proj$win, tperiod, m0=4.5, param0=param01, no.itr=1)
  res
 }
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
  \url{http://www.razi.ac.ir/jalilian/}
}
\keyword{spatial}
\keyword{print}
\keyword{models}

