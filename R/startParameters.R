#' Start Values for estimateCutoff
#'
#' A heuristic function to estimate start values for the ECOFFinder.
#'
#' @param DF data frame, containing columns "diam" and "obs" which specify the Zone Data
#' @param spMean method for estimating mean start parameter;  "mean" (default) or "peak1"
#'
#'
#' @return The function \code{startParameters} returns a numeric vector of estimated start parameters for mean, sd and n.
#'
#' @seealso \code{\link{estimateCutoff}}
#' @importFrom stats density weighted.mean

startParameters <- function(DF, spMean="mean"){
  # check input:
  stopifnot(spMean %in% c("mean", "peak1") & is.numeric(DF$diam) & is.numeric(DF$obs)) #numeric values are not null
  # n:
  n <- sum(DF$obs)
  # mean:
  mean <- numeric();   mx_dens <- numeric(); mx <- numeric()
  dens <- density(rep(DF$diam,times=DF$obs), bw=1)
  if(spMean=="mean" | sum(DF$obs)<30){mean <- weighted.mean(DF$diam, DF$obs)}
  else{
    mx_dens <- which(diff(diff(dens$y)>=0)<0)           #Find local Maxima
    mx <- tail(which(dens$y[mx_dens]>=max(dens$y[mx_dens]*0.125)), n=1)
    mean <- dens$x[mx_dens[mx]]
  }
  # sd:
  sd <- sqrt(Hmisc::wtd.var(DF$diam, weights=DF$obs))
  # combine and return:
  sp <-  c(mean = mean, sd =sd, n=n)
  names(sp)<- c("mean","sd", "n")
  return(sp)
}
