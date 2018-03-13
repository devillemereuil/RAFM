\name{shrews}
\Rdversion{1.1}
\alias{shrews}
\docType{data}
\title{
Genotype data from Finnish common shrews (Sorex araneus)
}
\description{
149 shrews from 5 subpopulations genotyped for 20 microsatellite markers
}
\usage{data(shrews)}
\format{
  A data frame with 149 observations on the following 21 variables.
  \describe{
    \item{\code{subpop}}{subpopulation code}
    \item{\code{X9}}{allelic copy 1, locus X9}
    \item{\code{X9.1}}{copy 2, locus X9}
    \item{\code{X14}}{copy 1, locus X14}
    \item{\code{X14.1}}{. . .}
    \item{\code{X33}}{. . .}
    \item{\code{X33.1}}{. . .}
    \item{\code{X45}}{. . .}
    \item{\code{X45.1}}{. . .}
    \item{\code{X62}}{. . .}
    \item{\code{X62.1}}{. . .}
    \item{\code{X67}}{. . .}
    \item{\code{X67.1}}{. . .}
    \item{\code{X68}}{. . .}
    \item{\code{X68.1}}{. . .}
    \item{\code{X69}}{. . .}
    \item{\code{X69.1}}{. . .}
    \item{\code{X92}}{. . .}
    \item{\code{X92.1}}{. . .}
    \item{\code{X97}}{. . .}
    \item{\code{X97.1}}{. . .}
  }
}
\details{
	Each row represents one individual.
	In this case, the genotypes are microsatellite fragment lengths.
	Missing data is coded as NA.
	}
\source{
Hanski I & Kuitunen J (1986): Shrews on small islands: epigenetic variation elucidates population stability. Holarctic Ecology 9: 193-204.
}

\examples{
data(shrews)
table(c(shrews[,2:3])) # allele counts of locus X9
}
\keyword{datasets}
