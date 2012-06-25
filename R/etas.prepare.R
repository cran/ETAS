
etas.prepare <- function(X, win, tperiod, m0, bwd, nnp, bwm, verbose)
{
  if (!is.numeric(tperiod) || length(tperiod) !=2)
    stop(paste(sQuote(tperiod), "must be a numeric vector of length 2."))
  tstart <- X$domain$xrange[1]
  tend <- X$domain$xrange[2]
  tstart2 <- tperiod[1]
  tlength <- tperiod[2]
  if (tstart2 < tstart || tstart2 > tend || tlength + tstart > tend)
    stop(paste(sQuote(tperiod), "is not inside the time domain [", tstart, ",", tend, "]."))
  tstart2 <- tstart2 - tstart

  if (!is.numeric(m0) || length(m0) != 1 || m0 <= 0 )
    stop(paste("the magnitude threshold", sQuote(m0), "must be a positive mumber"))

  verifyclass(win, "owin")
  switch(win$type, polygonal= {
    px <- win$bdry[[1]]$x
    py <- win$bdry[[1]]$y
   }, rectangle = {
    px <- c(win$xrange, win$xrange[2:1])
    py <- rep(win$yrange, each=2)
  })
  # repeat the first vertex
  np <- length(px) + 1
  px[np] <- px[1]
  py[np] <- py[1]
  rpoly <- cbind(px, py)

  verifyclass(X, "ppx")
  ctype <- X$ctype
  cnames <- names(X$data)
  tname <- cnames[ctype == "temporal"]
  sname <- cnames[ctype == "spatial"]
  mname <- cnames[ctype == "mark"]
  data <- as.list(X$data)
  xx <- data[[sname[1]]]         # x-coordinate of spatial coordinates
  yy <- data[[sname[2]]]         # y-coordinate of spatial coordinates
  tt <- data[[tname]]            # temporal coordinate: time
  mm <- data[[mname[1]]]         # mark values: magnitude

  ok <- (tt <= tlength + tstart) & (tt > tstart) & (mm >= m0)
  xx <- xx[ok]
  yy <- yy[ok]
  tt <- tt[ok] - tstart
  mm <- mm[ok] - m0
  flag <- as.integer(inside.owin(xx, yy, win))
  flag[tt < tstart2] <- -2
  revents <- cbind(tt, xx, yy, mm, flag, 0, 1, 0)

  if (is.null(bwd))
  {
    rbwd <- nndist.default(xx, yy, k=nnp)
    rbwd <- pmax(rbwd, bwm)
  }
  else
    rbwd <- bwd

  if (verbose) 
  cat("\tnumber of events: ", nrow(revents), 
      "\n\tnumber of target events: ", sum(flag == 1), 
      "\n\t[tstart, tend]: \t[", tstart , ", " , tend, "]", 
      "\n\ttstart2: ", tstart2, "\ttlength", tlength,
      "\n\tmagnitude threshold: ", m0, "\n")

  out <- list(revents=revents, rpoly=rpoly, rbwd=rbwd, 
              rtperiod=c(tstart2, tlength), sname=sname, tname=tname)
  return(out)
}

