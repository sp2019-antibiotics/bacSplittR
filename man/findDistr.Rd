% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findDistr.R
\name{findDistr}
\alias{findDistr}
\title{Fitting a Distribution}
\usage{
findDistr(df, startParams, fit, q)
}
\arguments{
\item{df}{data frame, containing columns "diam" and "obs" which specify the Zone Data}

\item{startParams}{numeric vector with start parameters for the nonlinear least squares algorithm}

\item{fit}{criterion for the fit; default: "n_abs"; options: "Pmean", "Psd", "Pn", "sigma"}

\item{q}{numeric value, giving the quantile of the distribution to be used as the ECOFF estimate}
}
\value{
A numeric vector with the estimated ECOFF, mean, standard deviation, n and the fit value.
}
\description{
Fitting a n-fold multiple normal distribution with the parameters mean, standard deviation (sd) and n to a given subset of Data.
}
\details{
The technique used is non-linear least squares, so the sum of the squared residuals is minimized.
}
\seealso{
\code{\link{estimateCutoff}}
}
