#' Estimate Epidemiological Cutoff Values
#'
#' The function estimateCutoff contains an algorithm designed to estimate an Epidemiological Cutoff Value of Zone Data pairings of bacteria and antibiotics.
#' The method was inspired by the paper "Statistical characterisation of bacterial wild-type MIC value distributions and the determination of epidemiological cut-off values".
#'
#'
#' @param obs numeric vector giving the absolute frequency of observations per diameter
#' @param diam numeric vector containing the diameter range (in mm), same length as obs, with default from 6:50
#' @param start method for estimating mean start parameter;  default: "mean", other option: "peak1"
#' @param fit criterium for the fit; default: "n_abs", other options: "Pmean", "Psd", "Pn", "sigma"
#' @param q numeric value for the quantile of the distribution; default: 0.0075
#' @param plot logical; if FALSE (default) no plot is drawn
#'
#' @return numeric value, estimate for the ECOFF
#'
#' @details
#'  Normal distributions are fitted to several subsets of the data, where (a part of) the most resistant data is removed.
#'  For this, Nonlinear Least Squares is used, which needs to be initzialized with starting values.
#'
#' There are two options to estimate the start value of the parameter mean of the fitted normal distribution:
#'\itemize{
#'  \item{"mean" = weighted mean}
#'  \item{"peak1" = the rightmost peak of the kernel density estimate, which surpasses a certain share (12.5%) of the mode.
#'  This can be especially useful, if the resistant observations also resemble a normal distribution. }
#' }
#'
#' To determine the subset which yields the best fitted normal distribution, one of the following criteria needs to be chosen:
#'\itemize{
#'\item n_abs: the absolute deviation of the estimated number of observations in the subset to N
#'\item Pmean: the p-value of the parameter mean in the NLS model
#'\item Psd: the p-value of the parameter sd in the NLS model
#'\item Pn: the p-value of the parameter n in the NLS model
#'\item sigma: a parameter of the nls-function, which is the square root of the estimated variance of the random error
#'}
#'
#' Here, the default value is "n_ab", which is also used in the original algorithm for MIC data.
#'
#' @examples
#' data("ZD", package = "EUCASTData") #load data
#' observations <- as.numeric(ZD[706,4:48])
#' estimateCutoff(observations, start = "peak1", fit = "n_abs", plot=TRUE)
#' example1 <- as.numeric(subset(ZD, Antimicrobial == "Ampicillin" & Bacterium == "Escherichia coli",
#'   grepl("^Z", colnames(ZD))))
#' estimateCutoff(example1, start = "peak1", fit = "n_abs", plot=TRUE)
#' example2 <- as.numeric(subset(ZD, Antimicrobial == "Piperacillin" & Bacterium == "Escherichia coli",
#'   grepl("^Z", colnames(ZD))))
#' estimateCutoff(example2, start = "peak1", fit = "n_abs", plot=TRUE)
#' example3 <- as.numeric(subset(ZD, Antimicrobial == "Mecillinam" & Bacterium == "Escherichia coli",
#'   grepl("^Z", colnames(ZD))))
#' estimateCutoff(example3, start = "peak1", fit = "n_abs", plot=TRUE)
#'
#' @references Turnidge, J., Kahlmeter, G., Kronvall, G. (2006) Statistical characterization of
#'   bacterial wild-type MIC value distributions and the determination of
#'   epidemiological cut-off values. Clin Microbial Infect 12: 418-425
#'   doi: 10.1111/j.1469-0691.2006.01377.x
#'
#' @export
#' @importFrom utils tail

estimateCutoff <- function(
  obs,                  #number of observations per diameter
  diam = 6:50,          #diameter range (in mm), same length as obs
  start = "mean",       #method for estimating mean start parameter: mean/mode/peak1
  fit = "n_abs",        #which criterion for fitting should be used: Pmean/Psd/Pk/n_abs/sigma/...
  q = 0.0075,            #which quantile of the distribution should be used for the estimate
  plot = FALSE          #should a plot be drawn
){
  # check input
  stopifnot(length(obs)==length(diam),
            start %in% c("mean", "peak1"),
            fit %in% c("Pmean", "Psd", "Pn", "n_abs", "sigma"),
            is.numeric(q), length(q) == 1, q>0 & q<1,
            is.logical(plot))
  # generate dataframe with all necessary data
  DF <- data.frame(obs=obs, diam=diam, lower=diam-0.5)
  DF$s <- rev(cumsum(rev(obs)))

  # calculate start parameters
  startP <- startParameters(DF, spMean = start)

  # find count of iterations
  peak1 <- startParameters(DF, spMean="peak1")["mean"]           #maximal cut off range is limited by peak1
  imax <-  ifelse(peak1>(DF$diam[2]),sum(DF$diam<(peak1)),2)     #number of observations that could be removed
  m <- matrix(NA, ncol=5, nrow=imax)

  # compute all possible ecoffs and their fitting parameter:
  for(i in 1:imax){
    m[i,] <- findDistr(DF[i:nrow(DF),], startParams=startP, fit =fit, q=q)
  }

  # plot:
  if(plot & any(!is.na(m[,1]))) print(plotDistr(m, startP,DF))

  #return estimated cutoff value from best fit
  return(ifelse(any(!is.na(m[,5])),m[,1][which.min(m[,5])],NA))
}

