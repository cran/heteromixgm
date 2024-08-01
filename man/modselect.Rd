\name{modselect}
\alias{modselect}

\title{modselect}

\description{
Model selection using the AIC, BIC and eBIC.
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
\item{selectmat}{Matrix containing the "optimal" l1 and l2 values for each
information criterion.}
\item{theta_aic}{Estimated precision matrices using the AIC for model selection.}
\item{theta_bic}{Estimated precision matrices using the BIC for model selection.}
\item{theta_ebic}{Estimated precision matrices using the EBIC for model selection.}
}

\references{
1. Hermes, S., van Heerwaarden, J., & Behrouzi, P. (2024). 
Copula graphical models for heterogeneous mixed data. 
Journal of Computational and Graphical Statistics, 1-15. \cr
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
