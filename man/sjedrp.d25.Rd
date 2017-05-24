\name{sjedrp.d25}
\alias{sjedrp.d25}
\title{Example data set for DRP (dmin 25).}
\description{
  Example data sets for checking the DRP code.  These data sets have
  been run on Rodieck's program for comparison.
}
\usage{data(sjedrp.d25)}


\format{Two columns of data; first column is x coordinate; second column
  is y coordinate.}

\source{These data sets are artificial; they were created with the dmin rule.}

\references{
  Rodieck RW (1991). The density recovery profile: A method for the
  analysis of points in the plane applicable to retinal studies. 
  Visual Neuroscience 6: 95-111.
}
\seealso{
  `sjedrp' for the main code,  and   `rodieckach' for the data to
  generate the DRP from Rodieck (1992).  See sjedrp.d10 for
  more documentation.
}

\examples{
data(sjedrp.d25)
bb <-c(0, 700, 0, 700)
res <- autodrp(sjedrp.d25$x, sjedrp.d25$y, nbins=20, r=10, a=bb)
plot(res)
stopifnot(abs(0.28 - res$p)< 0.03)
stopifnot(abs(23.0 - res$effrad)< 0.3)

}

\author{Stephen Eglen}
\keyword{datasets}

