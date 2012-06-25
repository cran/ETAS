\name{lambda}
\alias{lambda}
\title{Clustering Part of Conditional Intensity Function of the ETAS Model}
\description{
  A function to compute the clustering part of the conditional 
  intensity function of the ETAS model at specified time and location. 
}
\usage{
  lambda(t, x, y, param, X, m0)
}
\arguments{
  \item{t}{A numeric value. The time that the conditional intensity is 
            to be computed at.}
  \item{x}{A numeric value. The x-coordinate of the location that the 
           conditional intensity is to be computed at.}
  \item{y}{A numeric value. The y-coordinate of the location that the 
           conditional intensity is to be computed at.}
  \item{param}{Vector of model paramters.}
  \item{X}{Marked spatio-temporal point pattern (object of class
    \code{ppx}) representing time, sptial location, magnitude and
    ... of earthquakes.}
  \item{m0}{Minimum magnitude threshold.}
}
\value{
  A numeric value.
}
\details{
  For a given \eqn{t}{t}, \eqn{x}{x} and \eqn{y}{y}, this function 
  computes
    \deqn{
      \sum_{t_i < t} k(m_i)g(t - t_i)f(x - x_i, y - y_i|m_i).}{sum[t[i] < t] k(m[i]) g(t - t[i]) f(x - x[i], y - y[i]|m[i]).}
}
\references{
  Zhuang, J., Ogata, Y. and Vere-Jones, D. (2005). 
  Diagnostic analysis of space-time branching processes for earthquakes. 
  Lecture Note in Statistics: \emph{Case Studies in Spatial Point Process Models} 
  (Baddeley, A., Gregori, P., Mateu, J., Stoica, R. and Stoyan, D.), 
  Springer-Verlag, New York, \bold{185}, 276--292.

  Zhuang, J., Ogata, Y. and Vere-Jones, D. (2002). 
  Stochastic declustering of space-time earthquake occurrences. 
  \emph{Journal of the American Statistical Association},
  \bold{97}, 369--380.
}

\seealso{
  \code{\link{etas}}

  \code{\link{ppx}}
}
\examples{
  param00 <- c(0.592844590,
    0.204288231,
    0.022692883,
    1.495169224,
    1.109752319,
    0.001175925,
    1.860044210,
    1.041549634)

  win <- owin(c(41, 66), c(24, 42))
  proj <- long2flat(iran.quakes, win)
  lambda(50435, 40.12, 34.5, param00, proj$X, 4.5)
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
  \url{http://www.razi.ac.ir/jalilian/}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
