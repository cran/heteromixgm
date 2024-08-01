\name{plot_pcorgraph}
\alias{plot_pcorgraph}

\title{Plot partial correlation graphs}

\description{
Plots all \eqn{K} partial correlation graphs based on the \eqn{\Theta} selected 
using one of the information criteria.
}

\usage{
plot_pcorgraph(Theta, pos_clr, neg_clr, plot_layout, label_cex)
}

\arguments{
  \item{Theta}{List of \eqn{K} selected \eqn{\Theta}}
  \item{pos_clr}{Color, hexadecimal color allowed, representing the positive
  partial correlations in the plotted graphs.}
  \item{neg_clr}{Color, hexadecimal color allowed, representing the negative
  partial correlations in the plotted graphs.}
  \item{plot_layout}{Number of rows and columns for the plot layout.}
  \item{label_cex}{Size of the vertex labels in the plotted graphs.}
}

\value{
There is no return value. The function only shows plots in the graphics output device.
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
temp <- data_sim(network = "Random", n = 100, p = 20, K = 4, ncat = 6, rho = 0.25,
         gamma_o = 0.5, gamma_b = 0.1, gamma_p = 0.2, prob = 0.05)
X <- temp$z
l1 <- c(0.1)
l2 <- c(0,0.1)
gamma <- 0.5
ncores <- 1
est <- heteromixgm(X, "Approximate", l1, l2, ncores)
temp = modselect(est, X, l1, l2, gamma)
plot_pcorgraph(temp$theta_aic, "green", "red", c(2,2), 4.5)
}
}
