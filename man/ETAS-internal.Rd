\name{ETAS-internal}
\title{Internal ETAS functions}
\alias{etas.prepare}
\alias{rates}
\alias{lambdax}
\alias{etasfit}
\alias{decluster}
\description{
  Internal ETAS functions.
}
\usage{
etas.prepare(X, win, tperiod, m0, bwd, nnp, bwm, verbose)
rates(theta, revents, win, tperiod, bwd, dimyx=NULL, method="zhuang", plot.it=TRUE)
lambdax(rt, rx, ry, theta, revents)
etasfit(theta, revents, rpoly, tperiod, integ0, verbose)
decluster(theta, rbwd, revents, rpoly, tperiod)
}

\details{
  These are usually not to be called by the user.
}
\keyword{internal}
