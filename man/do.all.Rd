\name{do.all}
\Rdversion{1.1}
\alias{do.all}
\title{
Estimating coefficients of coancestry by admixture F-model
}
\description{
FST values and coancestry matrix by push of a button.
}
\usage{
do.all(dat, nMC, burnin, thinning, eps=7, priq=1, pria=c(1,2), prik=1, dom=FALSE, unobs=FALSE)
}

\arguments{
  \item{dat}{
data matrix, format as in dataset specimen
}
  \item{nMC}{
total number of MCMC iterations
}
  \item{burnin}{
no. of burnin iterations to be cast away, adjustment of proposal distributions is done during burnin
}
  \item{thinning}{
interval to sample from the MCMC, e.g. thinning=10 samples every 10th iteration
}
  \item{eps}{
truncation threshold of the Dirichlet distributions. eps=7 truncates at 1e-07, eps=8 at 1e-08, etc.
}

  \item{priq}{
specifies prior for ancestral allele frequencies, Dirichlet(priq) 
}

  \item{pria}{
specifies a prior for log(alpha): first component mean, second variance
}

  \item{prik}{
specifies a prior for kappa
} 

  \item{dom}{
are the data dominant (e.g. AFLP, default is FALSE)
}

  \item{unobs}{
is there a probability of unobserved alleles (e.g. micro-satellites, default is FALSE)
}

}

\details{ 
This function performs all functionality of RAFM in one go, i.e. it's a wrapper of AFM and the subsequent parameter tranformations gen.theta & gen.fst. The output is a list which contains samples from the posterior distributions of population-genetic parameters. The user needs to play around with AFM and/or different values of input parameters, if convergence problems occur.
}

\value{
  \item{theta}{Posterior of coancetry matrix theta. Dimension is (population x population x mc iteration).}
  \item{fst}{Posterior of Wright's fixation index. Each number is an mc iteration.} 
  \item{kap}{Posterior of lineage loadings kappa. Dimension is (population x lineage x mc iteration).}
  \item{alpha}{Posterior of drift parameters alpha. Dimension is (lineage x mc iteration).}
  }

\author{
Markku Karhunen, Uni. Helsinki
}

\examples{
data(specimen)
samp <- do.all(specimen, 100, 50, 2)
}
