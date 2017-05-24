\name{sjedrp.d10}
\alias{sjedrp.d10}
\title{Example data set for DRP (dmin 10).}
\description{
  Example data sets for checking the DRP code.  These data sets have
  been run on Rodieck's program for comparison.
}
\usage{data(sjedrp.d10)}


\format{Two columns of data; first column is x coordinate; second column
  is y coordinate.}

\source{These data sets are artificial; they were created with the dmin rule.}

\references{
  Rodieck RW (1991). The density recovery profile: A method for the
  analysis of points in the plane applicable to retinal studies. 
  Visual Neuroscience 6: 95-111.
}
\seealso{
  \link{autodrp} for the main code, and \link{rodieckach} for the data to
  generate the DRP from Rodieck (1992).  See other data sets
  such as \link{sjedrp.d20}, \link{sjedrp.d25} and \link{sjedrp.d30}.
}
\examples{
data(sjedrp.d10)
bb <- c(0, 700, 0, 700)
res <- autodrp(sjedrp.d10$x, sjedrp.d10$y, nbins=20, r=10, a=bb)
plot(res) 
stopifnot(abs(0.09 - res$p)< 0.03)
stopifnot(abs(12.9 - res$effrad)< 0.3)
}
\author{Stephen Eglen}
\keyword{datasets}
