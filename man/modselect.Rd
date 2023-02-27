\name{modselect}
\alias{modselect}

\title{modselect}

\description{
Model selection using the AIC and eBIC.
}

\usage{
modselect(est, X, l1, l2, gamma)
}

\arguments{
  \item{est}{Estimates of model obtained from cgmmd() function}
  \item{X}{A list of \eqn{K} \eqn{n_k} by \eqn{p} data matrices.}
  \item{l1}{Vector containing l1 penalty values.}
  \item{l2}{Vector containing l2 penalty values.}
  \item{gamma}{EBIC gamma parameter.}
}

\value{
\item{aic_idx}{Index of the estimated model corresponding to the optimal model
as per the AIC.}
\item{ebic_idx}{Index of the estimated model corresponding to the optimal model
as per the eBIC.}
\item{l1_aic}{Optimal l1 value as per the AIC.}
\item{l2_aic}{Optimal l2 value as per the AIC.}
\item{l1_ebic}{Optimal l1 value as per the eBIC.}
\item{l2_ebic}{Optimal l1 value as per the eBIC.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2022). Copula graphical
models for heterogeneous mixed data. arXiv preprint, arXiv:2210.13140. \cr
}

\author{
  Sjoerd Hermes, Joost van Heerwaarden and Pariya Behrouzi\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
\donttest{
X <- list(matrix(runif(25), 5, 5),matrix(runif(25), 5, 5),matrix(runif(25),
5, 5))
l1 <- c(0.4)
l2 <- c(0,0.1)
gamma <- 0.5
ncores <- 1
est <- heteromixgm(X, "Approximate", l1, l2, ncores)
modselect(est, X, l1, l2, gamma)
}
}
