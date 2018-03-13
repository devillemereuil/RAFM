\name{AFLP}
\Rdversion{1.1}
\alias{AFLP}
\docType{data}
\title{
A small simulated AFLP dataset.
}
\description{
Format your own AFLP data according to this sample dataset for use in AFM. 
}
\usage{data(AFLP)}
\format{
  A data frame with 100 observations on the following 21 variables.
  \describe{
    \item{\code{subpop}}{subpopulation code}
    \item{\code{loc1}}{AFLP band 1}
    \item{\code{loc2}}{AFLP band 2}
    \item{\code{loc3}}{AFLP band 3}
    \item{\code{loc4}}{. . .}
    \item{\code{loc5}}{. . .}
    \item{\code{loc6}}{. . .}
    \item{\code{loc7}}{. . .}
    \item{\code{loc8}}{. . .}
    \item{\code{loc9}}{. . .}
    \item{\code{loc10}}{. . .}
    \item{\code{loc11}}{. . .}
    \item{\code{loc12}}{. . .}
    \item{\code{loc13}}{. . .}
    \item{\code{loc14}}{. . .}
    \item{\code{loc15}}{. . .}
    \item{\code{loc16}}{. . .}
    \item{\code{loc17}}{. . .}
    \item{\code{loc18}}{. . .}
    \item{\code{loc19}}{. . .}
    \item{\code{loc20}}{. . .}
  }
}
\details{
Each row represents one individual. Column 1 contains the subpopulation code (can be any number or string). Columns 2+ contain the presence/absence information for all AFLP bands. Presences should be coded 1 and absence 0. Missing data must be coded by NA. Subpopulation codes must NOT be missing, and there has to be at least two subpopulations. (Insert a dummy row, if you want to estimate within-population theta for one subpopulation.) 
}

\examples{
data(AFLP)
head(AFLP) # how it looks
write.table(AFLP, row.names=FALSE, file="AFLP.txt") # for reference
}
\keyword{datasets}
