
lambda <- function(t, x, y, param, X, m0)
{
  if (!is.numeric(t) || !is.numeric(x) || !is.numeric(y))
    stop(paste("Arguments", sQuote(t), ",", sQuote(x), "and", sQuote(y), "must be numeric."))
  if (length(t) != length(x) || length(t) != length(y) || length(x) != length(y))
    stop(paste("Arguments", sQuote(t), ",", sQuote(x), "and", sQuote(y), "must be of the same length."))

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
  mm <- data[[mname[1]]] - m0    # mark values: magnitude
  revents <- cbind(tt, xx, yy, mm)
  storage.mode(revents) <- "double"
  theta <- sqrt(param)
  out <- numeric(length(t))
  for (i in 1:length(t))
  {
    if (t[i] < revents[1,1])
      out[i] <- 0
    else
    {
      out[i] <- .Call("lambdax", as.double(t[i]), as.double(x[i]), as.double(y[i]), 
                 as.double(theta), revents, PACKAGE="ETAS")[[1]]
    }
  }
  return(out)
}

