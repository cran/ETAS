
long2flat <- function(X, win)
{
  ymean <- centroid.owin(win)$y
  theta <- cos(ymean * pi / 180)
  A <- matrix(c(theta, 0, 0, 1), ncol=2, nrow=2)
  win <- affine(win, mat=A, vec=c(0,0), rescue=TRUE)
  xname <- names(X$data)[X$ctype == "spatial"][1]
  X$data[, xname] <- theta * as.list(X$data)[[xname]]
  X$domain$yrange <- theta * X$domain$yrange
  return(list(X = X, win = win, ymean=ymean))
}

