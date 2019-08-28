#' Check fit of nls estimate
#'
#'
#' @param M  nls object (or NA)
#' @param fit the method used to determine the best fit; options: "n_abs" (default), "Pmean", "Psd", "Pn", "sigma"
#' @param N the actual (true) number of observations in the subset (numeric value)
#'
#'
#' @return
#' Returns a numeric value describing the fit, this should be minimized to choose a subset
#'
#' @details
#' The parameter fit has the following options:
#'\itemize{
#'\item n_abs: the absolute deviation of the estimated number of observations in the subset to N
#'\item Pmean: the p-value of the parameter mean in the NLS model
#'\item Psd: the p-value of the parameter sd in the NLS model
#'\item Pn: the p-value of the parameter n in the NLS model
#'\item sigma: a parameter of the nls-function, which is the square root of the estimated variance of the random error
#'}
#'
#' @seealso \code{\link{estimateCutoff}}
#' @importFrom stats nls nls.control coefficients
#' @importFrom methods is

checkFit <- function(M, fit,N){ #M ...  nls object
  # check input:
  stopifnot(fit %in% c("n_abs", "Pmean", "Psd", "Pn", "sigma"), is.numeric(N), length(N) == 1, N>0)
  # choose measure corresponding to fit parameter:
  if (!is(M, "nls")) {
      result <- NA
  } else {
      result <- switch(fit,
                       Pmean = summary(M)$coefficients[1, 4],
                       Psd = summary(M)$coefficients[2, 4],
                       Pn = summary(M)$coefficients[3, 4],
                       n_abs = abs(coefficients(M)[3] - N),
                       sigma = summary(M)$sigma)
  }
  result
}
