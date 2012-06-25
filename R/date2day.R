
date2day <- function(dates, start, tz="", ...)
{
  start <- as.POSIXlt(start, tz=tz, ...)
  dates <- as.POSIXlt(dates, tz=tz, ...)
  out <- as.numeric(difftime(dates, start, units="days"))
  return(out)
}

