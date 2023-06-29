\name{heteromixgm}
\alias{heteromixgm}

\title{heteromixgm}

\description{
This function implements either the Gibbs or approximation method within the
Gaussian copula graphical model to estimate the conditional expectation for
the data that not follow Gaussianity assumption (e.g. ordinal, discrete,
continuous non-Gaussian, or mixed dataset).
}
\usage{
heteromixgm(X, method, lambda1, lambda2, ncores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A list containing \eqn{K} \eqn{n_k \times p} matrices (\eqn{K} is the number of groups, \eqn{n_k} is
  the sample size for group \eqn{k} and \eqn{p} is the number of variables}
  \item{method}{Choice between "Gibbs" and "Approximate" indicating which
  method to use.}
  \item{lambda1}{Vector containing values (in [0,1]) for the sparsity
  penalization of each \eqn{\Theta^k}.}
  \item{lambda2}{Vector containing values (in [0,1]) for the similarity
  penalization between the \eqn{\Theta^k}.}
  \item{ncores}{Number of cores to be used during parallel computing.}
}

\value{
\item{Z}{ New transformation of the data based on given or default
\code{Sigma}.}
\item{ES}{Expectation of covariance matrix( diagonal scaled to 1) of the
Gaussian copula graphical model.}
\item{Sigma}{The covariance matrix of the latent variable given the data.}
\item{Theta}{The inverse covariance matrix of the latent variable given the
data.}
\item{loglik}{Value of the Log likelihood under the estimated parameters.}
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
data(maize)
l1 <- c(0.4)
l2 <- c(0,0.1)
ncores <- 1
est <- heteromixgm(maize, "Approximate", l1, l2, ncores)
}
}

