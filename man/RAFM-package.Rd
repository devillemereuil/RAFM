\name{RAFM-package}
\Rdversion{1.2}
\alias{RAFM-package}
\alias{RAFM}
\docType{package}
\title{ Admixture F-model by Metropolis-Hastings
	}

\description{ This package involves codes needed for fitting the AFM from command line, conversion of AFM to coancestry and FST, and three example datasets. * * New feature in version 1.2: you can use function do.all to do all analyses! * *

	}
\details{
\tabular{ll}{
Package: \tab RAFM\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2013-10-24\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
In versions 1.2+, the primary mode to use this package is to run do.all with sufficient sample size.

If convergence problems occur, you should first increase the sample size, and then look into lower level, i.e. into AFM and gen.theta. If that does not help, you are free to play around with the source code (available at the software page, EGRU, University of Helsinki).

The code can be very slow for polymorphic data, such as shrews and sticklebacks in this package, so please be patient.  
}
\author{
Markku Karhunen, <markku.karhunen@helsinki.fi>
}
\references{
Karhunen & Ovaskainen (2012): Estimating population-level coancestry coefficients by an admixture F-model. Genetics 192: 609-617.
}
\examples{
data(specimen)
samp <- do.all(specimen, 100, 50, 2)
}
