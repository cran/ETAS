
check.temporal <- function(X)
{
   tname <- names(X$data)[X$ctype == "temporal"]
   tt <- as.list(X$data)[[tname]]
   for (i in 2:nrow(X$data))
      if ( tt[i-1] > tt[i])
         stop(paste("event", sQuote(i), "must be replaced with event", sQuote(i-1)))

  cat("OK! The temporal coordinate of the pattern is in ascending order.\n")
}

