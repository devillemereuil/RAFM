\name{heat.theta}
\alias{heat.theta}

\title{
Illustrating neutral genetic differentiation
}
\description{
This function illustrates the posterior of coancestry matrix theta as a heat map.
}
\usage{
heat.theta(thetapost, mean=TRUE, verbose=TRUE)
}

\arguments{
  \item{thetapost}{
posterior of theta
}
  \item{mean}{
mean=TRUE uses posterior means, mean=FALSE uses posterior medians
}
  \item{verbose}{
verbose=TRUE returns you summary statistics
}
}
\details{
This very simple function calculates the posterior mean/median of the coancestry matrix theta, and uses R's standard heatmap to illustrate it. Do not blame me for the dendogram or choice of colors. More sophisticated illustrations can be produced by manipulating the posterior manually.

Ugly averaging of off-diagonals takes care that the output matrices are positive definite, as they should for a sound theta.
}
\value{
If verbose=TRUE, you will get a list that has three components: posterior mean/median, and 95 percent posterior credible interval of the theta matrix.
}

\examples{
# Example commented away for obvious reasons
# data(specimen)
# samp <- do.all(specimen, 2000, 1000, 2) # <- obvious
# stats <- heat.theta(samp$theta, mean=FALSE, verbose=TRUE)
# names(stats)
# round(stats$median, 2) # and copy-paste to manuscript

# I do not claim responsibility for the dendogram!
}
