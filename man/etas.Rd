\name{etas}
\alias{etas}
\title{Fit the space-time ETAS model to data}
\description{
  A function to fit the space time epidemic type aftershock sequence 
  (ETAS) model to a catalog (point pattern) of earthquakes and perform
  a stochastic declustering method. 
}
\usage{
 etas(X, win, tperiod, m0, param0, bwd = NULL, nnp = 5, bwm = 0.05, 
      verbose = TRUE, plot.it = TRUE, no.itr = 11)
}
\arguments{
  \item{X}{Marked spatio-temporal point pattern (object of class
    \code{"ppx"}) representing time, spatial location, magnitude and
    ... of earthquakes.}
  \item{win}{Target geographical region of study (object of class
    \code{"owin"}). Events outside this region are treated as
    complementary events in the log-likelihood.}
  \item{tperiod}{Time period of quakes. A numeric vector of size 2: 
        (tstart, tlength). \code{tstart} is the starting point of 
         time target interval and \code{tlength} is the total length 
         of time period.}
  \item{m0}{Minimum magnitude threshold. A positive numeric value.}
  \item{param0}{Initial guess for model parameters. A numeric vector 
                of appropriate length (currently 8).}
  \item{bwd}{Optional. Bandwidths for smoothness and integration 
             on the geographical region \code{win}. A numeric vector 
             which has the length of the number of events.
             If not supplied, the following arguments \code{nnp} and 
             \code{bwm} determine bandwidths.}
  \item{nnp}{Number of nearest neighbors for bandwidth calculations. An integer.}
  \item{bwm}{Minimum bandwidth. A positive numeric value.}
  \item{verbose}{Logical flag indicating whether to print progress reports.}
  \item{plot.it}{Logical flag indicating whether plot probabilities of
    each event being a background even on a map.}
  \item{no.itr}{An integer indicating the number of iterations for 
               convergence of the iterative approach of declustering algorithm. See details.}
}
\value{
  A list with components
  \describe{
    \item{param:}{The ML estimates of model parameters.}
    \item{bk:}{An estimate of the \eqn{u(x, y)}{u(x, y)}.}
    \item{pb:}{The probabilities of being background event.}
    \item{opt:}{The results of optimization: the value of the log-likelihood
      function at the optimum point, its gradient at the optimum point and AIC of the model.}
    \item{rates:}{Pixel images of the estimated total intensity, 
          background intensity, clustering intensity and conditional intensity.}
  }
}
\details{
  Ogata (1988) introduced the epidemic type aftershock sequence 
  (ETAS) model based on Gutenberg-Richter law and modified Omori law. 
  In its space-time representation (Ogata, 1998), the ETAS model is a 
  temporal marked point process model, and a special case of marked 
  Hawks process, with conditional intensity function
    \deqn{
      \lambda(t, x, y | H_t) = \mu(x,y) + \sum_{t_i < t} k(m_i)g(t - t_i)f(x - x_i, y - y_i|m_i)
    }{
      lambda(t, x, y | H_t ) = mu(x, y) 
       + sum[t[i] < t] k(m[i]) g(t - t[i]) f(x - x[i], y - y[i]|m[i])
    }
  where
  \describe{
      \item{\eqn{H_t}{H_t}: }{is the observational history up to time t, but not 
            including t; that is
            \deqn{H_t=\{(t_i, x_i, y_i, m_i): t_i < t\}}{H_t = {(t[i], x[i], y[i], m[i]): t[i] < t}}
           }
      \item{\eqn{\mu(x,y)}{mu(x, y)}: }{is the background intensity. Currently it 
            is assumed to take the semi-parametric form
            \deqn{\mu(x,y)=\mu u(x,y)}{mu(x, y) = mu u(x, y)}
            where \eqn{\mu}{mu} is an unknown constant and \eqn{u(x,y)}{u(x, y)} 
            is an unknown function.
           }
      \item{\eqn{k(m)}{k(m)}: }{is the expected number of events triggered 
            from an event of magnitude \eqn{m}{m} given by
           \deqn{k(m) = A\exp(\alpha(m - m_0))}{k(m[i]) =  A exp(alpha(m - m0))}
           }
      \item{\eqn{g(t)}{g(t)}: }{is the p.d.f of the occurrence 
            times of the triggered events, taking the form
           \deqn{g(t) =  \frac{p-1}{c}(1 + \frac{t}{c})^{-p}}{g(t) = ((p - 1)/c)(1 + t/c)^(-p)}
           }
      \item{\eqn{f(x,y|m)}{f(x, y|m)}: }{is the p.d.f of 
            the locations of the triggered events, considered to be
            either the long tail inverse power density
           \deqn{ f(x, y|m) = \frac{q-1}{\pi \sigma(m))}
             (1 + \frac{x^2 + y^2}{\sigma(m)})^{-q} }{ f(x, y|m) = (q - 1)/(pi sigma(m)) 
                 (1 + (x^2 + y^2)/(sigma(m)))^(-q)}
           or the light tail Gaussian density (currently not implemented)
           \deqn{ f(x,y|m)= \frac{1}{2\pi \sigma(m)}\exp(-\frac{x^2 + y^2}{2\sigma(m)}) }{f(x,y|m)= exp(-(x^2 + y^2)/(2 sigma(m)))/(2 pi sigma(m))}
           with
           \deqn{ \sigma(m) = D\exp(\gamma(m - m_0)) }{ sigma(m) = D exp( gamma (m - m0) )} 
       } 
  }

  The ETAS models classify seismicity into two components, background 
  seismicity \eqn{\mu(x, y)}{mu(x,y)} and clustering seismicity 
  \eqn{\lambda(t, x, y|H_t) - \mu(x, y)}{lambda(t, x, y|H_t) - mu(x, y)}, where 
  each earthquake event, whether it is a background event or generated by 
  another event, produces its own offspring according to the branching rules
  controlled by \eqn{k(m)}{k(m)}, \eqn{g(m)}{g(m)} and \eqn{f(x, y|m)}{f(x, y|m)}.

  Background seismicity rate \eqn{u(x, y)}{u(x, y)} and the model parameters 
  \deqn{\theta=(\mu, A, c, \alpha, p, D, q, \gamma)}{theta = (mu, A, c, alpha, p, D, q, gamma)}
  are estimated simultaneously using an iterative approach proposed in Zhuang et al. (2002). 
  First, for an initial \eqn{u_0(x, y)}{u0(x, y)}, the parameter vector 
  \eqn{\theta}{theta}  is estimated by maximizing the log-likelihood function 
  \deqn{l(\theta)=\sum_{i} \lambda(t_i, x_i, y_i|H_{t_i}) - \int \lambda(t, x, y|H_t) dx dy dt.}{l(theta) = sum_[i] lambda(t[i], x[i], y[i]|H_t) 
         - int lambda(t, x, y|H_t) dx dy dt.} 
  Then the procedure calculates the probability of being a background 
  event for each event in the catalog by
  \deqn{ \phi_i = \frac{\mu(x_i, y_i)}{\lambda(t_i, x_i, y_i|H_{t_i})}. }{mu(x[i], y[i])/lambda(t[i], x[i], y[i]|H_t[i]).}
  Using these probabilities and kernel smoothing method with Gaussian kernel 
  and appropriate choice of bandwidth (determined by \code{bwd} or \code{nnp} 
  and \code{bwm} arguments), the background rate \eqn{u_0(x, y)}{u0(x, y)} 
  is updated. These steps are repeated enough times 
  (determined by \code{no.itr} argument) such that the 
  results converge.

  The procedure is a \code{C} port of the original \code{Fortran} code 
  by Ogata, Zhunang and their colleagues.
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

  Ogata, Y. (1998). 
  Space-time point-process models for earthquake occurrences. 
  \emph{Annals of the Institute of Statistical Mathematics}, 
  \bold{50}, 379--402.

  Ogata, Y. (1988). 
  Statistical models for earthquake occurrences and residual analysis for point processes.
  \emph{Journal of American Statistical Association}, 
  \bold{83}, 9--27.
}
\seealso{

  \code{\link{ppx}} and \code{\link{owin}} for constructing data.

  \code{\link{etasfit}}.
}
\examples{
  # fitting the ETAS model to an Iranian catalog
  \dontrun{
  plot(iran.quakes)
  }
  # specifiying the geographical region
  win <- owin(c(41, 66), c(24, 42))
  \dontrun{
  plot(iran.quakes$data[, 2:3])
  plot(win, add=TRUE)
  }
  # projecting log-lat coordinates into flat map coordinates
  proj <- long2flat(iran.quakes, win)
  # specifying time period
  tperiod <- c(25000, 40329)
  # initial parameters values
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
  }



  # fitting the ETAS model to a Japanese catalog
  \dontrun{
  plot(jap.quakes)
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
  tperiod <- c(10000, 23376)
  # initial parameters values
  param00 <- c(0.592844590,
   0.204288231,
   0.022692883,
   1.495169224,
   1.109752319,
   0.001175925,
   1.860044210,
   1.041549634)
  \dontrun{
  res <- etas(proj$X, proj$win, tperiod, m0=4.5, param0=param00, no.itr=11)
  }
}
\author{Abdollah Jalilian
  \email{jalilian@razi.ac.ir}
  \url{http://www.razi.ac.ir/jalilian/}
}
\keyword{spatial}
\keyword{math}
\keyword{earthquake modeling}
