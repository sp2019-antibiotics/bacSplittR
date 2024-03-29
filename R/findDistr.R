#' Fitting a Distribution
#'
#' Fitting a n-fold multiple normal distribution with the parameters mean, standard deviation (sd) and n to a given subset of Data.
#'
#' @param df data frame, containing columns "diam" and "obs" which specify the Zone Data
#' @param fit  criterion for the fit; default: "n_abs"; options: "Pmean", "Psd", "Pn", "sigma"
#' @param startParams numeric vector with start parameters for the nonlinear least squares algorithm
#' @param q numeric value, giving the quantile of the distribution to be used as the ECOFF estimate
#'
#' @return A numeric vector with the estimated ECOFF, mean, standard deviation, n and the fit value.
#'
#'
#' @details The technique used is non-linear least squares, so the sum of the squared residuals is minimized.
#'
#' @seealso \code{\link{estimateCutoff}}
#' @importFrom stats coefficients qnorm nls nls.control
#' @importFrom methods is

findDistr <- function(df, startParams, fit, q){
  # check input:
  stopifnot(
   is.numeric(df$diam), is.numeric(df$s),
   is.numeric(startParams), length(startParams)==3,
   fit %in% c("n_abs", "Pmean", "Psd", "Pn", "sigma"),
   is.numeric(q), length(q) == 1, q>0 & q<1
  )
  # fit n-fold multiple of normal survival function:
  M <- tryCatch(
      nls(s~n * pnorm(lower, mean=mean, sd=sd, lower.tail = FALSE), data=df, start=startParams,
          lower=c(min(df$diam),0,0), algorithm="port", control=nls.control(maxiter=10000)),
      error=function(e){NA})

  # estimate cutoff:
  e <- if (is(M, "nls")) qnorm(q,coefficients(M)[1], coefficients(M)[2]) else NA
  # if(!is.na(e)& e<(min(df$diam)-0.5)) e <- NA

  # estimated coefficients:
  est <- if (is(M, "nls")) coefficients(M)[1:3] else c(NA, NA, NA)

  # calculate fit measurement
  g <-  checkFit(M,fit,max(df$s))

  # combine and return:
  res <- c(e,est,g)
  names(res) <- c("estimate","estMean","estSD","estN","fit")
  return(res)
}
