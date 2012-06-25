
rates <- function(theta, revents, win, tperiod, bwd, dimyx=NULL, method="zhuang", plot.it=TRUE)
{
   xx <- revents[, 2]
   yy <- revents[, 3]
   mm <- revents[, 4]
   flag <- revents[, 5]
   bk <- revents[, 6]
   pb <- revents[, 7]
   lam <- revents[, 8]
   switch(method, spatstat={
     Xg <- ppp(xx, yy, window=owin(range(xx), range(yy)))
     total.im <- density(Xg, weights=rep(1/diff(tperiod), Xg$n), dimyx)[win]
     bkgd.im <-  smooth.ppp(Xg %mark% bk, weights=bk/diff(tperiod))[win]
     clust.im <- smooth.ppp(Xg %mark% (lam - bk), weights=rep(1/diff(tperiod), Xg$n))[win]
     lambd.im <- smooth.ppp(Xg %mark% lam)[win]
   }, zhuang={
     if (is.null(dimyx))
       dimyx <- c(128, 128)
     if (!is.numeric(dimyx) || length(dimyx) != 2)
       stop(paste(sQuote(dimyx), "must be a numeric vector of length 2."))
     gr <-gridcenters(bounding.box(win),  dimyx[2], dimyx[1])
     gx <- gr$x
     gy <- gr$y
     out1 <- out2 <- out3 <- out4 <- numeric(length(gx)) 
     for (i in 1:length(gx))
       {
          r2 <- (xx - gx[i])^2 + (yy - gy[i])^2
          s1 <- exp(-r2/(2 * bwd^2)) / (2 * pi * bwd^2)
          s2 <- pb *  s1
          s1 <- sum(s1)/diff(tperiod)
          s2 <- sum(s2)/diff(tperiod)
          out1[i] <- s1
          out2[i] <- s2
          out3[i] <- s1 - s2
          lamb <- theta[1] * s2 + lambdax(tperiod[2], gx[i], gy[i], theta, revents)
          out4[i] <- lamb
      }
     total.im <- as.im(list(x=unique(gx), y=unique(gy), z=matrix(out1, nrow=dimyx[1], ncol=dimyx[2])))
     bkgd.im <- as.im(list(x=unique(gx), y=unique(gy), z=matrix(out2, nrow=dimyx[1], ncol=dimyx[2])))
     clust.im <- as.im(list(x=unique(gx), y=unique(gy), z=matrix(out3,, nrow=dimyx[1], ncol=dimyx[2])))
     lambd.im <- as.im(list(x=unique(gx), y=unique(gy), z=matrix(out4, nrow=dimyx[1], ncol=dimyx[2])))
   })
   if (plot.it)
   {
     dev.new()
     par(mfrow=c(2, 2), mar=c(1.5, 1.5, 2, 3))
     plot(total.im, main="total spatial seismicity rate", axes=TRUE)
     plot(bkgd.im, main="background seismicity rate", axes=TRUE)
     plot(clust.im, main="clustering rate", axes=TRUE)
     plot(lambd.im, main="conditional intensity", axes=TRUE)
   }
   return(list(total.im, bkgd.im, clust.im, lambd.im))
}

