% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tune_delta.R
\name{tune_delta}
\alias{tune_delta}
\title{Tune delta for momentLS estimators}
\usage{
tune_delta(x, nSplits = 1, c_M_const = 0, fcn.summarize = "mean")
}
\arguments{
\item{x}{input chain}

\item{nSplits}{number of splits}

\item{c_M_const}{constant in \eqn{c_M = c_M_const * \sqrt{\log M}} in the mhat criterion in 
\eqn{\hat{m} = \min\{t\in 2\mathbb{N}; \hat{\rho}_M(t+2) \le c_M \sqrt{\log M/M}\}}}

\item{fcn.summarize}{if fcn.summarize =="mean", the average of deltas from each split will be returned. If fcn.summarize =="median", median of deltas from each split will be returned.}
}
\value{
A numeric value (delta) and a data frame (mhat_all)
\itemize{
  \item delta: estimated delta (mean or median from \eqn{\hat{delta}_M} from each split.)
  \item mhat_all: estimated \eqn{\hat{m}} and \eqn{\hat{delta}_M} from each split.
}
}
\description{
Tune based on a modification of an adaptive bandwidth selection method by Politis <https://doi.org/10.1080/10485250310001604659>. 
See also our paper <https://arxiv.org/abs/2207.12705>.
}
