
etas <- function(X, win, tperiod, m0, param0, bwd = NULL, nnp = 5, bwm = 0.05, 
                 verbose = TRUE, plot.it = TRUE, no.itr = 11)
{
  prep <- etas.prepare(X, win, tperiod, m0, bwd, nnp, bwm, verbose)
  revents <- prep$revents
  rpoly <- prep$rpoly
  rbwd <- prep$rbwd
  rtperiod <- prep$rtperiod
  sname <- prep$sname
  tname <- prep$tname
  N <- nrow(revents)
  k <- ncol(revents)

  if (!is.numeric(param0) || length(param0) != 8 || any(param0 < 0))
    stop(paste(param0), "must be a numeric vector of length 8 with positive components.")


  param1 <- param0
  names(param1) <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  thetar <- matrix(, nrow=no.itr, ncol=8)

  for (itr in 1:no.itr)
  {
    bkg <- decluster(param1, rbwd, revents, rpoly, tperiod)
    revents <- bkg$revents
    integ0 <- bkg$integ0
    bk <- revents[,6]
    pb <- revents[,7]
    if (verbose)
    {
      cat("iteration: ", itr, "\n")
      cat("======================================================\n")
      cat("background seismicity rate:\n")
      print(summary(bk))
      cat("probability of being a background event:\n")
      print(summary(pb))
      cat("integral of background seismicity rate: ", integ0, "\n")
      cat("======================================================\n")
    }
    if (plot.it)
    {
      par(mfrow=c(1, 2), mar=c(4, 4, 3, 1))
      cols <- ifelse(pb < 0.5, "red", "blue")
      plot(revents[,2], revents[,3], cex = 0.05 + 2.5 * revents[,4]/m0, col=cols,
           main=paste("iteration: ", itr), xlab=sname[1], ylab=sname[2])
      plot(win, add=TRUE)
      plot(revents[,1], pb, xlab=tname, ylab="probability of being a background event")
      rates0 <- rates(param1, revents, win, tperiod, rbwd, plot.it=plot.it)
    }
    opt <- etasfit(param1, revents, rpoly, tperiod, integ0, verbose)
    thetar[itr,] <- opt$estimate
    param1 <- thetar[itr,]
    if (verbose)
    {
      cat("======================================================\n")
      cat("MLE:\n")
      print(param1)
      cat("======================================================\n")
    }
  }

  rates0 <- rates(param1, revents, win, tperiod, rbwd, plot.it=plot.it)
  names(param1) <- c("mu", "A", "c", "alpha", "p", "D", "q", "gamma")
  out <- list(param = param1, bk=bk, pb=pb, opt=opt, rates=rates0)
  class(out) <- "etas"
  return(out)
}


print.etas <- function (x, ...) 
{
  cat("ETAS model: fitted using iterative stochastic declustering method\n")
  cat("ML estimates of model parameters:\n")
  print(x$param)
  cat("log-likelihood: ", x$opt$loglik, "\nAIC: ", x$opt$aic, "\n")
}

#setMethod("print", signature="etas", print.etas,
#          where = topenv(parent.frame()),
#          valueClass = NULL, sealed = FALSE)

